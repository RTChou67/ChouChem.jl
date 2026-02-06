using Printf
using Test

# --- 1. Define the ChouChem Results (Hardcoded from previous run for stability) ---
# Format: (Molecule, Basis, Charge, Method, ChouChem_E)
ChouChem_Results = [
    ("HF", "STO-3G", "0", "RHF", -98.4147504766),
    ("HF", "STO-3G", "0", "RMP2", -98.4731521696),
    ("HF", "STO-3G", "0", "RCI", -98.5193602375),
    ("HF", "STO-3G", "0", "RCCSD", -98.5121725785),
    ("HF", "STO-3G", "1", "UHF", -98.1167714125),
    ("HF", "STO-3G", "1", "UCI", -98.1416672333),
    ("HF", "STO-3G", "1", "UCCSD", -98.1331607149),
    ("HF", "6-31G", "0", "RHF", -99.8550549975),
    ("HF", "6-31G", "0", "RMP2", -100.0142172859),
    ("HF", "6-31G", "0", "RCI", -100.0142910642),
    ("HF", "6-31G", "0", "RCCSD", -100.0016716547),
    ("HF", "6-31G", "1", "UHF", -99.4226454172),
    ("HF", "6-31G", "1", "UCI", -99.5203612370),
    ("HF", "6-31G", "1", "UCCSD", -99.5121950991),
    ("H2O", "STO-3G", "0", "RHF", -74.8492756164),
    ("H2O", "STO-3G", "0", "RMP2", -74.9219031622),
    ("H2O", "STO-3G", "0", "RCI", -74.9522844807),
    ("H2O", "STO-3G", "0", "RCCSD", -74.9523733132),
    ("H2O", "STO-3G", "1", "UHF", -74.5963826117),
    ("H2O", "STO-3G", "1", "UCI", -74.6657235576),
    ("H2O", "STO-3G", "1", "UCCSD", -74.6581125626),
    ("H2O", "6-31G", "0", "RHF", -75.8727852965),
    ("H2O", "6-31G", "0", "RMP2", -76.0331924707),
    ("H2O", "6-31G", "0", "RCI", -76.0327560547),
    ("H2O", "6-31G", "0", "RCCSD", -76.0293597201),
    ("H2O", "6-31G", "1", "UHF", -75.4417696906),
    ("H2O", "6-31G", "1", "UCI", -75.5597838978),
    ("H2O", "6-31G", "1", "UCCSD", -75.5511527163)
]

# --- 2. Function to parse Gaussian Logs ---
function parse_gaussian_log(filename)
    if !isfile(filename)
        return NaN
    end
    
    content = read(filename, String)
    
    # 1. Look for CCSD first (most specific)
    # "E(Corr)= -0.213...     E(TOT)= -76.029..."
    # Or "E(CCSD)= ..."
    # Gaussian CCSD output: "E(CCSD) =  -0.75029359720D+02"
    m_ccsd = match(r"E\(CCSD\)\s*=\s*([-0-9.D+eE]+)", content)
    if m_ccsd !== nothing
        return parse(Float64, replace(m_ccsd.captures[1], "D" => "e"))
    end

    # 2. Look for CISD
    # "E(CISD) =  -0.76032756055D+02"
    m_cisd = match(r"E\(CISD\)\s*=\s*([-0-9.D+eE]+)", content)
    if m_cisd !== nothing
        return parse(Float64, replace(m_cisd.captures[1], "D" => "e"))
    end

    # 3. Look for MP2
    # "EUMP2 =  -0.76033192471D+02"
    m_mp2 = match(r"EUMP2\s*=\s*([-0-9.D+eE]+)", content)
    if m_mp2 !== nothing
        return parse(Float64, replace(m_mp2.captures[1], "D" => "e"))
    end

    # 4. Look for SCF (RHF/UHF)
    # "SCF Done:  E(RHF) =  -75.8727852965"
    m_scf = match(r"SCF Done:\s+E\([RU]HF\)\s*=\s*([-0-9.]+)", content)
    if m_scf !== nothing
        return parse(Float64, m_scf.captures[1])
    end

    return NaN
end

# --- 3. Comparison Loop ---

println("\n" * "="^110)
println(@sprintf("%-8s | %-8s | %-3s | %-8s | %-16s | %-16s | %-16s", "Mol", "Basis", "Q", "Method", "ChouChem", "Gaussian", "Diff (Hartree)"))
println("="^110)

correct_count = 0
total_count = 0

for (mol, basis, q, method, chou_e) in ChouChem_Results
    filename = joinpath("test_g16", "$(mol)_$(basis)_Q$(q)_$(method).log")
    
    # Check if filename needs adjustment?
    # Test script generated: "$(mol_name)_$(basis)_Q$(charge)_$(method_name).gjf"
    # Matches exactly.
    
    g16_e = parse_gaussian_log(filename)
    
    diff = abs(chou_e - g16_e)
    
    # Tolerance: 1e-5 usually acceptable for comparisons, but we expect tight agreement for RHF/MP2.
    # CI/CCSD might differ slightly due to frozen core? No, STO-3G/6-31G usually all electron?
    # Gaussian defaults to Frozen Core for post-HF!
    # ChouChem RunRMPn / RunUCI / RunRCCSD might be Full or Frozen?
    # RunRMPn documentation says "Restricted MPn". Code uses all orbitals usually unless specified.
    # RunUCI uses active space. If MaxExcitation=2, it uses all singles/doubles.
    # Gaussian default is "FC" (Frozen Core). 
    # IF DIFF IS LARGE (> 1e-3), it's likely FC vs Full.
    # Gaussian keyword for Full is "Full". Our inputs were just "#p Method/Basis".
    
    status = diff < 1e-5 ? "OK" : "DIFF"
    if isnan(g16_e) status = "MISSING" end
    
    if status == "OK"
        global correct_count += 1
    end
    global total_count += 1

    println(@sprintf("%-8s | %-8s | %-3s | %-8s | %16.10f | %16.10f | %16.10e %s", mol, basis, q, method, chou_e, g16_e, diff, status))
end
println("="^110)
println("Match Rate: $correct_count / $total_count")
