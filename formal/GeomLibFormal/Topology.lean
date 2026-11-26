/-!
  GeomLibFormal — combinatorial topology core.

  Core-Lean only. A triangulated closed surface is modelled as a dart set
  with an edge involution `alpha` and a face permutation `phi` whose cycles
  have length 3. The vertex permutation is `sigma = phi ∘ alpha`.

  Two facts are proved: sigma is invertible (with inverse `alpha ∘ phi^2`),
  and the boundary-of-boundary of a triangle cancels over Z2 (checked over a
  finite vertex domain with `decide`).
-/

set_option maxRecDepth 1000000
set_option maxHeartbeats 8000000

namespace GeomLibFormal.Topology

/-- Vertex permutation induced by a combinatorial map. -/
def sigma {n : Nat} (alpha phi : Fin n → Fin n) (d : Fin n) : Fin n :=
  phi (alpha d)

def Involutive {α : Sort u} (f : α → α) : Prop := ∀ x, f (f x) = x

/-- If alpha is an involution and phi^3 = id, alpha ∘ phi^2 is a left and
 * right inverse of sigma = phi ∘ alpha. -/
theorem sigma_inverse {n : Nat} (alpha phi : Fin n → Fin n)
    (ha : Involutive alpha) (hphi : ∀ d, phi (phi (phi d)) = d) :
    (∀ d, alpha (phi (phi (sigma alpha phi d))) = d) ∧
    (∀ d, sigma alpha phi (alpha (phi (phi d))) = d) := by
  constructor
  · intro d
    dsimp [sigma]
    rw [show phi (phi (phi (alpha d))) = alpha d from hphi (alpha d)]
    exact ha d
  · intro d
    dsimp [sigma]
    rw [ha (phi (phi d))]
    exact hphi d

/-- Boundary of a triangle: its three directed edges. -/
def triangleEdges (a b c : Nat) : List (Nat × Nat) := [(a,b), (b,c), (c,a)]

/-- Over Z2 the vertex boundary of an edge ignores direction. -/
def edgeVertices (e : Nat × Nat) : List Nat := [e.1, e.2]

/-- Every vertex of a triangle appears an even number of times in the
 * boundary of its boundary. Checked exhaustively over a finite domain. -/
theorem boundary_of_boundary_even : ∀ a b c : Fin 4,
    let vs := (triangleEdges a.val b.val c.val).flatMap edgeVertices
    ((List.filter (· == a.val) vs).length % 2 == 0) &&
    ((List.filter (· == b.val) vs).length % 2 == 0) &&
    ((List.filter (· == c.val) vs).length % 2 == 0) := by
  decide

end GeomLibFormal.Topology
