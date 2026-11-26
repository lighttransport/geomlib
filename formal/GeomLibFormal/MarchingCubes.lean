/-!
  GeomLibFormal — Marching cubes case combinatorics.

  Core-Lean only (no mathlib dependency): every theorem is checked by the
  Lean kernel with `decide` over the finite type `Fin 256`.

  These facts mirror the invariants required by geomlib's marching-cubes
  implementation: each cube face has an even number of crossing edges (so
  neighbouring cubes share a consistent number of intersections), the total
  crossing count has the same parity as the number of inside corners, and
  complementing all signs preserves the crossing-edge set (complement
  symmetry).
-/

set_option maxRecDepth 1000000
set_option maxHeartbeats 8000000

namespace GeomLibFormal.MarchingCubes

/-- First endpoint of cube edge e (corner bits: x=1, y=2, z=4). -/
def aNat (e : Fin 12) : Nat :=
  match e.val with
  | 0 => 0 | 1 => 2 | 2 => 4 | 3 => 6
  | 4 => 0 | 5 => 1 | 6 => 4 | 7 => 5
  | 8 => 0 | 9 => 1 | 10 => 2 | 11 => 3
  | _ => 0

/-- Second endpoint of cube edge e. -/
def bNat (e : Fin 12) : Nat :=
  match e.val with
  | 0 => 1 | 1 => 3 | 2 => 5 | 3 => 7
  | 4 => 2 | 5 => 3 | 6 => 6 | 7 => 7
  | 8 => 4 | 9 => 5 | 10 => 6 | 11 => 7
  | _ => 0

/-- Bit k of an 8-bit sign mask. -/
def sign (mask : Fin 256) (k : Nat) : Bool := (mask.val / (2 ^ k)) % 2 == 1

/-- An edge crosses the isosurface iff its endpoints have different signs. -/
def crosses (mask : Fin 256) (e : Fin 12) : Bool :=
  sign mask (aNat e) != sign mask (bNat e)

def crossingEdges (mask : Fin 256) : List (Fin 12) :=
  (List.finRange 12).filter (crosses mask)

/-- Coordinate of corner k along axis 0=x, 1=y, 2=z. -/
def coord (k : Nat) (axis : Nat) : Nat :=
  match axis with
  | 0 => k % 2
  | 1 => (k / 2) % 2
  | _ => (k / 4) % 2

/-- Face f: 0/1 = x=0/1, 2/3 = y=0/1, 4/5 = z=0/1. -/
def onFace (k : Nat) (f : Fin 6) : Bool :=
  coord k (f.val / 2) == f.val % 2

def crossingOnFace (mask : Fin 256) (f : Fin 6) : List (Fin 12) :=
  (List.finRange 12).filter (fun e =>
    crosses mask e && onFace (aNat e) f && onFace (bNat e) f)

/-- Complement all eight corner signs. -/
def complement (mask : Fin 256) : Fin 256 :=
  ⟨(255 - mask.val) % 256, Nat.mod_lt _ (by decide)⟩

/-- Every cube edge joins two corners differing in exactly one bit. -/
theorem endpoints_differ_in_one_bit : ∀ e : Fin 12,
    let b0 := ((aNat e / 1) % 2 == (bNat e / 1) % 2)
    let b1 := ((aNat e / 2) % 2 == (bNat e / 2) % 2)
    let b2 := ((aNat e / 4) % 2 == (bNat e / 4) % 2)
    aNat e != bNat e && ((b0 && b1 && !b2) || (b0 && !b1 && b2) || (!b0 && b1 && b2)) := by
  decide

def insideCount (mask : Fin 256) : Nat :=
  ((List.finRange 8).filter (fun k : Fin 8 => sign mask k.val)).length

/-- Handshaking on the cube graph: crossing-count parity equals the parity of
 * the number of inside corners (each inside corner has odd degree 3). -/
theorem crossing_parity : ∀ mask : Fin 256,
    (crossingEdges mask).length % 2 = insideCount mask % 2 := by
  decide

theorem even_crossings_per_face : ∀ (mask : Fin 256) (f : Fin 6),
    (crossingOnFace mask f).length % 2 = 0 := by
  decide

theorem uniform_has_no_crossings :
    crossingEdges ⟨0, by decide⟩ = [] &&
    crossingEdges ⟨255, by decide⟩ = [] := by
  decide

theorem complement_symmetry : ∀ mask : Fin 256,
    crossingEdges mask = crossingEdges (complement mask) := by
  decide

/-- The crossing-edge count cannot exceed the 12 cube edges. -/
theorem crossing_count_bounded : ∀ mask : Fin 256,
    (crossingEdges mask).length <= 12 := by
  decide

end GeomLibFormal.MarchingCubes
