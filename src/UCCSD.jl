export RunUCCSD

function RunUCCSD(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int; MaxIter = 128)
    println("--- Starting UCCSD Calculation ---")
    println("Running UHF...")
    UHF_Res = ChouChem.RunUHF(MolInAng, Charge, Multiplicity, MaxIter = MaxIter)
    
    return spin_orbital_ccsd(UHF_Res)
end