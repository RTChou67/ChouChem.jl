export RunRCCSD

function RunRCCSD(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int)
    println("--- Starting RCCSD Calculation ---")
    println("Running RHF...")
    RHF_Res = ChouChem.RunRHF(MolInAng, Charge, Multiplicity)
    
    println("Converting RHF to Spin-Orbitals...")
    UHF_Res = ChouChem.RHF2UHF(RHF_Res)
    
    return spin_orbital_ccsd(UHF_Res)
end