import Lake
open Lake DSL

package GeomLibFormal where
  leanVersion := "v4.36.0"

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.36.0"

@[default_target]
lean_lib GeomLibFormal where
  roots := #[`GeomLibFormal]
