/-!
  GeomLibFormal — exact predicate specifications.

  Core-Lean only. The predicates are specified as integer determinants: the C
  implementation computes the same signs via floating-point expansions.
  The equivalence of the homogeneous determinant and the cross-product form is
  checked exhaustively over a bounded integer domain with `decide`.
-/

set_option maxRecDepth 1000000
set_option maxHeartbeats 8000000

namespace GeomLibFormal.Predicates

/-- Determinant of a 3x3 integer matrix (row-major). -/
def det3 (m : Fin 3 → Fin 3 → Int) : Int :=
  m 0 0 * (m 1 1 * m 2 2 - m 1 2 * m 2 1)
  - m 0 1 * (m 1 0 * m 2 2 - m 1 2 * m 2 0)
  + m 0 2 * (m 1 0 * m 2 1 - m 1 1 * m 2 0)

/-- orient2d as the determinant of the homogeneous 3x3 matrix. -/
def orient2dDet (ax ay bx byc cx cy : Int) : Int :=
  det3 (fun i j =>
    match i.val, j.val with
    | 0, 0 => ax | 0, 1 => ay | 0, _ => 1
    | 1, 0 => bx | 1, 1 => byc | 1, _ => 1
    | _, 0 => cx | _, 1 => cy | _, _ => 1)

/-- The same orientation as the cross product of edge vectors. -/
def orient2dCross (ax ay bx byc cx cy : Int) : Int :=
  (bx - ax) * (cy - ay) - (byc - ay) * (cx - ax)

/-- Exhaustive finite-domain check of the two formulations. -/
theorem orient2d_det_agrees_finite :
    ∀ ax ay bx byc cx cy : Fin 4,
      orient2dDet ax.val ay.val bx.val byc.val cx.val cy.val =
      orient2dCross ax.val ay.val bx.val byc.val cx.val cy.val := by
  decide

def subRow (i : Fin 3) : Fin 4 := ⟨i.val + 1, Nat.succ_lt_succ i.isLt⟩

-- Removing column 0: minor column j maps to j+1 (same as subRow).
-- Removing column 1: j=0 -> 0, otherwise j+1.
-- Removing column 2: j<2 -> j, otherwise 3.
-- Removing column 3: j -> j.
def subColSkipSecond (j : Fin 3) : Fin 4 :=
  if j.val = 0 then ⟨0, by decide⟩ else ⟨j.val + 1, Nat.succ_lt_succ j.isLt⟩

def subColSkipThird (j : Fin 3) : Fin 4 :=
  if h : j.val < 2 then ⟨j.val, Nat.lt_trans h (by decide : 2 < 4)⟩
  else ⟨3, by decide⟩

def subColSkipFourth (j : Fin 3) : Fin 4 :=
  ⟨j.val, Nat.lt_trans j.isLt (by decide : 3 < 4)⟩

/-- Explicit incircle 4x4 determinant by Laplace expansion along the first row.
 * Rows are (x, y, x^2+y^2, 1); positive means inside the circumcircle of a
 * counter-clockwise triangle a,b,c. -/
def incircle (ax ay bx byc cx cy dx dy : Int) : Int :=
  let r : Fin 4 → Fin 4 → Int := fun i j =>
    match i.val with
    | 0 => match j.val with | 0 => ax | 1 => ay | 2 => ax*ax + ay*ay | _ => 1
    | 1 => match j.val with | 0 => bx | 1 => byc | 2 => bx*bx + byc*byc | _ => 1
    | 2 => match j.val with | 0 => cx | 1 => cy | 2 => cx*cx + cy*cy | _ => 1
    | _ => match j.val with | 0 => dx | 1 => dy | 2 => dx*dx + dy*dy | _ => 1
  r 0 0 * det3 (fun i j => r (subRow i) (subRow j))
  - r 0 1 * det3 (fun i j => r (subRow i) (subColSkipSecond j))
  + r 0 2 * det3 (fun i j => r (subRow i) (subColSkipThird j))
  - r 0 3 * det3 (fun i j => r (subRow i) (subColSkipFourth j))

end GeomLibFormal.Predicates
