/-!
  Generate include/geomlib/generated/mc_table.h from the verified marching
  cubes definitions in GeomLibFormal.MarchingCubes.

    lean formal/gen_mc_table.lean > include/geomlib/generated/mc_table.h

  The definitions are duplicated here so the generator runs with plain `lean`
  and no lake/mathlib setup.
-/

set_option maxRecDepth 1000000
set_option maxHeartbeats 8000000

def aNat (e : Fin 12) : Nat :=
  match e.val with
  | 0 => 0 | 1 => 2 | 2 => 4 | 3 => 6
  | 4 => 0 | 5 => 1 | 6 => 4 | 7 => 5
  | 8 => 0 | 9 => 1 | 10 => 2 | 11 => 3
  | _ => 0

def bNat (e : Fin 12) : Nat :=
  match e.val with
  | 0 => 1 | 1 => 3 | 2 => 5 | 3 => 7
  | 4 => 2 | 5 => 3 | 6 => 6 | 7 => 7
  | 8 => 4 | 9 => 5 | 10 => 6 | 11 => 7
  | _ => 0

def sign (mask : Fin 256) (k : Nat) : Bool := (mask.val / (2 ^ k)) % 2 == 1

def crosses (mask : Fin 256) (e : Fin 12) : Bool :=
  sign mask (aNat e) != sign mask (bNat e)

def crossingEdges (mask : Fin 256) : List (Fin 12) :=
  (List.finRange 12).filter (crosses mask)

def edgeMask (mask : Fin 256) : Nat :=
  (crossingEdges mask).foldl (fun acc e => acc + 2 ^ e.val) 0

#eval IO.println "/* Generated from formal/gen_mc_table.lean (Lean-verified case set)."
#eval IO.println " * Do not edit by hand; regenerate with `lean formal/gen_mc_table.lean`. */"
#eval IO.println "#ifndef GEOMLIB_MC_TABLE_H"
#eval IO.println "#define GEOMLIB_MC_TABLE_H"
#eval IO.println ""
#eval IO.println "#include <stdint.h>"
#eval IO.println ""
#eval IO.println "static const uint16_t geomlib_mc_edge_table[256] = {"
#eval for m in List.finRange 256 do
  IO.println ("    " ++ toString (edgeMask m) ++ ",")
#eval IO.println "};"
#eval IO.println ""
#eval IO.println "#endif /* GEOMLIB_MC_TABLE_H */"
