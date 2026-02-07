export CCSDResults

struct CCSDResults
	ENum::Int
	nocc::Int
	nvir::Int
	t1::Matrix{Float64}
	t2::Array{Float64, 4}
	Ecc::Float64
	Etot::Float64
end
