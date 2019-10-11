using NetworkMetaAnalysis
using Test

@testset "NetworkMetaAnalysis.jl" begin
    # Write your own tests here.
end






exit()



@time using NetworkMetaAnalysis

@time using PaddedMatrices, Random, StructuredMatrices
using StructuredMatrices: RaggedMatrix

emax(a, em, ed50) = em * a / (a + ed50)

function generate_mbnma_data(
    Ntrials, minarms, maxarms, total_treats = (4 + maxarms - minarms),
    em = 2 .+ 5 .* randn(total_treats), ed50 = randexp(total_treats) .+ randexp(total_treats) .+ randexp(total_treats) .+ randexp(total_treats)
)
    arms = minarms:maxarms
    treat = collect(1:total_treats)
    treatments = Int[]
    coloffset = 0
    coloffsets = Vector{Int}(undef, Ntrials)
    collengths = Vector{Int}(undef, Ntrials)
    for n in 1:Ntrials
        coloffsets[n] = coloffset
        cl = rand(arms)
        coloffset += cl
        collengths[n] = cl
        append!(treatments, sort!(@view(shuffle!(treat)[1:cl])))
    end
    responses = Vector{Float64}(undef, length(treatments))
    doses = similar(responses)
    i = 1
    for n in 1:Ntrials
        mun = 0.5randn()
        for _ in 1:collengths[n]
            t = treatments[i]
            dose = ed50[t] * exp(randn())#(randexp() + randexp()) * 5
            doses[i] = dose
            responses[i] = emax(dose, em[t], ed50[t]) + mun + 0.2randn()
            i += 1
        end
    end
    responses, RaggedMatrix(treatments,coloffsets,collengths,maxarms), doses, em, ed50
end
resp,treat,dose,em,ed50 = generate_mbnma_data(80, 2, 8, 13);

function study_v(study_counts)
    studies = similar(study_counts, sum(study_counts))
    ind = 1
    for i ∈ eachindex(study_counts)
        c = study_counts[i]
        for _ ∈ 1:c
            studies[ind] = i
            ind += 1
        end
    end
    studies
end


studies = study_v(treat.column_lengths); studies'
treat.data'
extrema(treat.data)
# include("/home/c285497/.julia/dev/NetworkMetaAnalysis/src/NetworkMetaAnalysis.jl")

gn, (d2,r2) = NetworkMetaAnalysis.GatherNetwork(studies, treat.data, dose, resp);

gn.α.offsets
gn.α.masks[1] |> bitstring
gn.α.masks[2] |> bitstring
gn.α.nreps
extrema(gn.α.offsets)

gn.δ.offsets
gn.δ.masks[1] |> bitstring
gn.δ.masks[2] |> bitstring
gn.δ.nreps'


