const LoadedBasisSets = Dict{String, Dict{Int, Vector{CGTF}}}()
const PROJECT_DIR = joinpath(@__DIR__, "..")
const BASIS_SET_DIR = joinpath(PROJECT_DIR, "Basis")
const BasisNameToFile = include(joinpath(BASIS_SET_DIR, "BasisList.jl"))

export generate_basis_list

function get_basis_set(name::String)
	if haskey(LoadedBasisSets, name)
		return LoadedBasisSets[name]
	end
	if !haskey(BasisNameToFile, name)
		error("Cannot find basis set '$name' in BasisList.jl.")
	end
	filename = BasisNameToFile[name]
	filepath = joinpath(BASIS_SET_DIR, filename)
	if !isfile(filepath)
		error("Cannot find basis set file: '$filepath'")
	end
	try
		basis_data = include(filepath)
		LoadedBasisSets[name] = basis_data
		return basis_data
	catch e
		error("Failed to load or parse basis set file '$filepath': $e")
	end
end



function generate_basis_list(molecule::Vector{Atom})
	BasisList = Basis[]
	println("--- Start Generating Basis Functions ---")
	for atom in molecule
		println("Processing atom '$(atom.symbol)', basis set '$(atom.basis_set)'...")
		basis_set_data = get_basis_set(atom.basis_set)
		if isnothing(basis_set_data)
			error("Cannot load basis set '$(atom.basis_set)' for atom '$(atom.symbol)', calculation aborted.")
		end
		if !haskey(basis_set_data, atom.Z)
			error("Cannot find data for atomic number $(atom.Z) ('$(atom.symbol)') in basis set '$(atom.basis_set)', calculation aborted.")
		end
		element_data = basis_set_data[atom.Z]
		for cgtf in element_data
			new_basis_function = Basis(
				cgtf.Type,
				cgtf.GTFs,
				atom.position,
			)
			push!(BasisList, new_basis_function)
		end
		println("  -> Added $(length(element_data)) basis functions for '$(atom.symbol)'.")
	end
	println("--- Basis function generation complete ---
")
	return BasisList
end