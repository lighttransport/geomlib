# AGENTS.md

Guidance for AI coding agents working in this repository.

For the project layout, build commands, testing instructions, architecture, and
coding conventions, see [CLAUDE.md](CLAUDE.md).

## Git Push Policy (mandatory pre-push checklist)

Pushing rewrites shared state. Before **any** `git push` (regular, force, or
`--force-with-lease`), an agent must complete every step below. Skipping a step
is a defect: once a large binary or leaked credential lands on `origin`, it is
in the public Git history even if a later commit removes it.

### 1. Audit the exact commits about to leave the machine

Always inspect the exact range you are pushing. For a regular push, use the
configured upstream; for a force or lease push, use the actual remote tip being
replaced:

```bash
# Regular push with a configured upstream.
RANGE="$(git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}')..HEAD"

# Force/lease push: use the remote tip explicitly instead.
# RANGE="origin/main..HEAD"

git rev-list --count "$RANGE"
git log --oneline "$RANGE"
```

Run all four checks below against **every** commit in that range, not just
`HEAD`. Any finding blocks the push until it is resolved and the complete audit
has been run again.

#### Check 1 — Credentials and sensitive data

Never push a commit that contains, or has ever contained at any point in its
history within the push range:

- API keys, bearer tokens, OAuth client secrets, AWS / GCP / Azure access keys,
  SSH private keys, PGP private keys, `.netrc`, `.env`, or `.npmrc` files with
  `_authToken`.
- Email/password pairs, JWTs, or Slack/Discord/GitHub webhook URLs.
- Internal hostnames, VPN configuration, or customer asset paths that are not
  already public.
- Personal user paths (`/home/<someone>/…`, `C:\Users\<someone>\…`) baked into
  source comments or commit messages.

Run the heuristic from the repository root against the exact push range:

```bash
RANGE="${RANGE:-$(git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}')..HEAD}"
git diff "$RANGE" -- . ':!**/*.md' ':!**/*.txt' \
  | grep -nEi 'AKIA[0-9A-Z]{16}|ASIA[0-9A-Z]{16}|AIza[0-9A-Za-z_-]{35}|(api[_-]?key|secret|token|password|passwd|bearer|private[_-]?key)[[:space:]]*[:=]' \
  || echo "No credential-shaped strings found in $RANGE"
```

This is a heuristic, not a guarantee. If a commit changes CI, package upload,
dependency download, or external-service integration code, inspect that diff
manually as well.

##### Automated secret scanners (preferred when available)

When `gitleaks` and/or `trufflehog` are installed, run them over the **same push
range**. They complement the grep heuristic but do not replace manual review.
Any finding blocks the push; treat the credential as compromised and rotate it
before scrubbing history. If a tool is missing, record that it was skipped.

```bash
RANGE="${RANGE:-$(git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}')..HEAD}"
BASE="$(git rev-parse "${RANGE%%..*}")"
REPO="$(git rev-parse --show-toplevel)"
BRANCH="$(git branch --show-current)"
# TruffleHog's go-git backend cannot read a worktree's .git file. Point it at
# the main worktree, whose .git is a directory, and select the branch.
MAIN_WT="$(git worktree list --porcelain | awk '/^worktree /{print $2; exit}')"

if command -v gitleaks >/dev/null 2>&1; then
  gitleaks detect --source "$REPO" --redact --log-opts="$RANGE" \
    && echo "gitleaks: clean over $RANGE" \
    || { echo "ABORT: gitleaks flagged secrets in $RANGE"; exit 1; }
else
  echo "gitleaks not installed; relying on grep and manual review."
fi

TRUFFLEHOG="${TRUFFLEHOG:-$HOME/go/bin/trufflehog}"
[ -x "$TRUFFLEHOG" ] || TRUFFLEHOG="$(command -v trufflehog || true)"
if [ -n "$TRUFFLEHOG" ] && [ -x "$TRUFFLEHOG" ]; then
  "$TRUFFLEHOG" git "file://$MAIN_WT" --since-commit="$BASE" \
    --branch="$BRANCH" --fail \
    && echo "trufflehog: clean over $RANGE" \
    || { echo "ABORT: trufflehog flagged secrets in $RANGE"; exit 1; }
else
  echo "trufflehog not installed; relying on grep and manual review."
fi
```

##### Personal paths and private asset names

Scan commit messages as well as added source lines. Private model/scene names
can leak through a performance note even when the diff itself is harmless.
Public history should describe the change generically, not name a customer or
proprietary asset.

```bash
RANGE="${RANGE:-$(git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}')..HEAD}"

# Personal paths in messages and source must be empty.
{ git log "$RANGE" --format='%B';
  git diff "$RANGE" -- 'CMakeLists.txt' 'src/*' 'include/*' 'tests/*' \
    'web/*' 'scripts/*' 'formal/*' 'cmake/*'; } \
  | grep -inE '/home/[a-z]|/Users/[A-Za-z]|C:\\Users\\' \
  && { echo "ABORT: personal path in messages or source"; exit 1; } \
  || echo "clean: no personal paths"

# Eyeball every asset-filename token in messages and added source lines.
{ git log "$RANGE" --format='%B';
  git diff "$RANGE" -- 'CMakeLists.txt' 'src/*' 'include/*' 'tests/*' \
    'web/*' 'scripts/*' 'formal/*' 'cmake/*' | grep '^+'; } \
  | grep -oiE '[A-Za-z0-9_-]{3,}\.(usdz|usda|usdc|png|jpg|jpeg|exr|gltf|glb|fbx|obj)' \
  | sort -u
```

Generic names and deliberately public fixtures are fine. Strip any private
asset name with an interactive rebase or `git filter-repo --replace-message`,
then re-run every Check 1 scan because rewriting changes commit hashes.

#### Check 2 — No build artifacts

Build outputs do not belong in Git. They inflate the repository, cause merge
conflicts on every rebuild, and often embed absolute paths from the developer's
machine. Reject paths such as:

- `build/`, `build-*/`, `web/build/`, `web/demo/dist/`, `node_modules/`.
- `*.a`, `*.o`, `*.obj`, `*.so`, `*.dylib`, `*.dll`, `*.lib`, `*.pdb`, `*.exe`,
  `*.out`.
- `*.wasm` and generated `web/demo/public/*.js` / `*.wasm` outputs.
- `*.ninja_deps`, `*.ninja_log`, `build.ninja`, `CMakeCache.txt`,
  `CMakeFiles/`, `CTestTestfile.cmake`, `cmake_install.cmake`, and
  `compile_commands.json`.
- `__pycache__/`, `*.pyc`, `*.egg-info/`, and other language build/cache state.
- `.claude/`, `.codex/`, `.commandcode/`, and other AI/tool scratch directories.

Run from the repository root:

```bash
RANGE="${RANGE:-$(git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}')..HEAD}"
git diff --name-only "$RANGE" \
  | grep -Ei '(^|/)(build[^/]*|CMakeFiles|__pycache__|node_modules|dist)(/|$)|\.(a|o|obj|so|dylib|dll|lib|pdb|pyc|wasm|ninja_deps|ninja_log)$|(^|/)(CMakeCache\.txt|CTestTestfile\.cmake|build\.ninja|compile_commands\.json|cmake_install\.cmake)$|^web/demo/public/.*\.(js|wasm)$' \
  && { echo "ABORT: build artifacts in $RANGE"; exit 1; } \
  || echo "No obvious build artifacts in $RANGE"
```

If an artifact appears, remove it from the affected commits and add an
appropriate `.gitignore` rule. Do not push it merely because a later commit
deletes it.

#### Check 3 — No unintended binary or asset data

GeomLib's tracked tree is source, tests, documentation, and small generated
source tables. Binary art, captured payloads, downloaded design documents, and
other non-regenerable data belong in external storage or Git LFS, not in the
main repository history.

Reject, unless the user explicitly approves a deliberate small fixture:

- `*.blend`, `*.fbx`, `*.glb`, `*.gltf`, `*.abc`, `*.mb`, `*.ma`.
- `*.usd`, `*.usda`, `*.usdc`, `*.usdz` captures outside an intentional fixture
  location.
- `*.png`, `*.jpg`, `*.jpeg`, `*.exr`, `*.hdr`, `*.tif`, or `*.tiff` over about
  256 KiB.
- Captured `.pdf` documents, archives, and files under capture/scratch paths.
- Anything that cannot be trivially regenerated from source.

Inspect the range's sizes and binary diffs:

```bash
RANGE="${RANGE:-$(git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}')..HEAD}"

# Binary blobs appear as "-" in the added/deleted columns.
git diff --numstat "$RANGE" \
  | awk '$1 == "-" || $2 == "-" { print "ABORT: binary diff:", $0; found=1 } END { exit found }' \
  && echo "No binary diffs in $RANGE" \
  || { echo "ABORT: binary diff in $RANGE"; exit 1; }

# Inspect every large change; text line counts are only a heuristic.
git diff --stat "$RANGE"

git diff --name-only "$RANGE" \
  | grep -iE '\.(blend|fbx|glb|gltf|abc|mb|ma|exr|hdr|tif|tiff|pdf)$|/captures?/|(^|/)(\.claude|\.codex|\.commandcode)/' \
  && { echo "ABORT: unexpected binary/asset paths in $RANGE"; exit 1; } \
  || echo "No flagged binary/asset paths in $RANGE"
```

If a binary genuinely must ship, use Git LFS or external storage and ask the
user before adding it.

#### Check 4 — User permission (always, no exceptions)

After Checks 1–3 pass, **stop and ask the user before pushing**, even if the
user previously said "push when done" or has approved many pushes in the same
session. Authorization is single-shot and scoped to the immediate request. The
required order is **audit, summarize, ask, then push**.

The summary must include:

- Branch and remote.
- Whether the push is fast-forward or a force push; for force pushes, list the
  previously public commit hashes that would be orphaned.
- Commit count and a one-line summary per commit.
- Results from Checks 1–3, including any scanner that was unavailable.
- Relevant test results and any known limitations.

Do not treat "auto mode" or a prior `--yes`-style flag as push authorization.
When the user approves a force push, prefer `--force-with-lease` over `--force`
so the push is rejected if the remote moved since the last fetch.

### 2. After an authorized push

- After a force push, report the previous remote tip so teammates can recover
  their work.
- Open or update the PR if one exists; GeomLib's default branch is `main`.
- Do not delete branches unless the user explicitly asks.

### 3. If a check fails

- **Credential leak:** stop, surface it, treat the credential as compromised,
  rotate first, scrub history second, and re-run the full checklist.
- **Build artifact or unintended asset:** remove it from every affected commit
  with an interactive rebase or `git filter-repo --invert-paths`, update
  `.gitignore`, and re-run the full checklist.
- **Local development note:** remove it from tracking with `git rm --cached`,
  add an ignore rule, and preserve the local file. `git rm --cached` does not
  erase earlier commits; if the note must never be published, rewrite the
  commits before pushing and re-run the full checklist.
- **Uncertain asset:** ask whether it belongs in Git LFS, external storage, or
  should be removed.

The checklist is non-negotiable. A single large binary permanently bloats every
future clone; a leaked key is a security incident regardless of how quickly it
is rotated.
