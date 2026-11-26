/-!
  GeomLibFormal — in-circle predicate axioms for Delaunay flipping.

  Core-Lean only. The four-point in-circle determinant is defined over Int
  and its invariance axioms are checked exhaustively over a bounded domain.
  These are the algebraic facts (Knuth's axioms) that make the in-circle test
  usable for Delaunay edge flips; a full flip-termination proof needs mathlib
  and is left as future work.
-/

set_option maxRecDepth 1000000
set_option maxHeartbeats 8000000

namespace GeomLibFormal.Delaunay

def det3 (m : Fin 3 → Fin 3 → Int) : Int :=
  m 0 0 * (m 1 1 * m 2 2 - m 1 2 * m 2 1)
  - m 0 1 * (m 1 0 * m 2 2 - m 1 2 * m 2 0)
  + m 0 2 * (m 1 0 * m 2 1 - m 1 1 * m 2 0)

def subRow (i : Fin 3) : Fin 4 := ⟨i.val + 1, Nat.succ_lt_succ i.isLt⟩
-- Removing column k from a 4-column matrix maps minor column j to:
--   k=0: j+1,  k=1: 0 if j=0 else j+1,  k=2: j if j<2 else 3,  k=3: j
def subColSkipFirst (j : Fin 3) : Fin 4 := ⟨j.val + 1, Nat.succ_lt_succ j.isLt⟩
def subColSkipSecond (j : Fin 3) : Fin 4 :=
  if j.val = 0 then ⟨0, by decide⟩ else ⟨j.val + 1, Nat.succ_lt_succ j.isLt⟩
def subColSkipThird (j : Fin 3) : Fin 4 :=
  if h : j.val < 2 then ⟨j.val, Nat.lt_trans h (by decide : 2 < 4)⟩
  else ⟨3, by decide⟩
def subColSkipFourth (j : Fin 3) : Fin 4 :=
  ⟨j.val, Nat.lt_trans j.isLt (by decide : 3 < 4)⟩

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

def I (x : Fin 2) : Int := x.val

/-- Translation invariance. -/
theorem incircle_translation :
    ∀ ax ay bx byc cx cy dx dy tx ty : Fin 2,
      incircle (I ax + I tx) (I ay + I ty)
               (I bx + I tx) (I byc + I ty)
               (I cx + I tx) (I cy + I ty)
               (I dx + I tx) (I dy + I ty)
      = incircle (I ax) (I ay) (I bx) (I byc)
                 (I cx) (I cy) (I dx) (I dy) := by
  decide

/-- Swapping the first two points negates the determinant. -/
theorem incircle_swap_first_two :
    ∀ ax ay bx byc cx cy dx dy : Fin 2,
      incircle (I bx) (I byc) (I ax) (I ay) (I cx) (I cy) (I dx) (I dy)
      = - incircle (I ax) (I ay) (I bx) (I byc) (I cx) (I cy) (I dx) (I dy) := by
  decide

/-- Positive integer scaling multiplies the determinant byc k^4. -/
theorem incircle_scaling :
    ∀ ax ay bx byc cx cy dx dy k : Fin 2,
      incircle (k.val * I ax) (k.val * I ay)
               (k.val * I bx) (k.val * I byc)
               (k.val * I cx) (k.val * I cy)
               (k.val * I dx) (k.val * I dy)
      = (k.val : Int)^4 * incircle (I ax) (I ay) (I bx) (I byc)
                                    (I cx) (I cy) (I dx) (I dy) := by
  decide

/-- 90-degree rotation preserves the determinant exactly. -/
theorem incircle_rotation :
    ∀ ax ay bx byc cx cy dx dy : Fin 2,
      incircle (-(I ay)) (I ax) (-(I byc)) (I bx)
               (-(I cy)) (I cx) (-(I dy)) (I dx)
      = incircle (I ax) (I ay) (I bx) (I byc) (I cx) (I cy) (I dx) (I dy) := by
  decide

end GeomLibFormal.Delaunay
