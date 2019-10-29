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
# include("~/.julia/dev/NetworkMetaAnalysis/src/NetworkMetaAnalysis.jl")

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







MetaAnalysis((fe1,fe2),gn,(fep1,fep2),(identity,exp))




using BenchmarkTools, SLEEFPirates
function uexp(v1, v2, v3, v4)
    v1 = SLEEFPirates.exp(v1)
    v2 = SLEEFPirates.exp(v2)
    v3 = SLEEFPirates.exp(v3)
    v4 = SLEEFPirates.exp(v4)
    (v1, v2, v3, v4)
end
function lexp(v1, v2, v3, v4)
    v = (v1, v2, v3, v4)
    @inbounds for i ∈ eachindex(v)
        v = Base.setindex(v, SLEEFPirates.exp(v[i]), i)
    end
    v
end
function vaexp(v...)
    @inbounds ntuple(i -> SLEEFPirates.exp(v[i]), length(v))
end
v1 = ntuple(Val(4)) do _ Core.VecElement(randn()) end
v2 = ntuple(Val(4)) do _ Core.VecElement(randn()) end
v3 = ntuple(Val(4)) do _ Core.VecElement(randn()) end
v4 = ntuple(Val(4)) do _ Core.VecElement(randn()) end

uexp(v1, v2, v3, v4)
lexp(v1, v2, v3, v4)
vaexp(v1, v2, v3, v4)

@benchmark uexp($v1, $v2, $v3, $v4)
@benchmark lexp($v1, $v2, $v3, $v4)
@benchmark vaexp($v1, $v2, $v3, $v4)


@code_native uexp(v1, v2, v3, v4)
@code_native lexp(v1, v2, v3, v4)



using SIMDPirates, VectorizationBase
const W64, Wshift64 = VectorizationBase.pick_vector_width_shift(Float64)


function expand_test!(expanded, mul, a, nreps)
    ind = 1
    @inbounds for i ∈ eachindex(a)
        aᵢ = a[i]
        for j ∈ 1:nreps[i]
            expanded[ind] = aᵢ
            ind += 1
        end
    end
    t = vbroadcast(Vec{W64,Float64}, 0.0)
    ptre = pointer(expanded); ptrm = pointer(mul)
    for i ∈ 0:(length(expanded) >> Wshift64)-1
        t = vfmadd(vload(Vec{W64,Float64}, ptre), vload(Vec{W64,Float64}, ptrm), t)
        ptre += 8W64
        ptrm += 8W64
    end
    t
end
function gather_test(mul, a, offs)
    t = vbroadcast(Vec{W64,Float64}, 0.0)
    ptra = pointer(a); ptrm = pointer(mul); ptri = pointer(offs)
    for i ∈ 0:(length(mul)>>Wshift64)-1
        vm = vload(Vec{W64,Float64}, ptrm)
        va = gather(vadd(ptra, vload(Vec{W64,Int}, ptri)))
        t = vfmadd(vm, va, t)
        ptrm += 8W64
        ptri += sizeof(Int)*W64
    end
    t
end

function calc_offs_from_nreps(nreps)
    offs = Vector{Int}(undef, sum(nreps))
    ind = 1
    for i ∈ eachindex(nreps)
        for _ ∈ 1:nreps[i]
            offs[ind] = (i - 1) * sizeof(Int)
            ind += 1
        end
    end
    offs
end

nreps = rand(2:4, 50);
a = randn(50);
offs = calc_offs_from_nreps(nreps);
mul = randn(length(offs));
expanded = similar(mul);

expand_test!(expanded, mul, a, nreps)
gather_test(mul, a, offs)

using BenchmarkTools

@benchmark expand_test!($expanded, $mul, $a, $nreps)
@benchmark gather_test($mul, $a, $offs)

@code_native gather_test(mul, a, offs)



