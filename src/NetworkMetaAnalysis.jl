module NetworkMetaAnalysis

using SIMDPirates, StackPointers, ReverseDiffExpressionsBase, PaddedMatrices

export NetworkMetaAnalysis, RaggedNetwork, GatherNetwork

abstract type AbstractNetwork end

"""
Represents a network as a ragged matrix.
"""
struct RaggedNetwork <: AbstractNetwork

end
"""
Represents a network as a series of offsets.
"""
struct GatherNetwork <: AbstractNetwork
    arms_and_studies::Matrix{Int}
    masks::Vector{UInt8}
    narms::Vector{Int}
end
function GatherNetwork(studies, arms)
    @boundscheck length(studies) == length(arms) || throw("There should be one study ID per arm.")
    if !issorted(studies)
        studies_sp = sortperm(studies)
        studies = studies[studies_sp]
        arms = arms[studies_sp]
    end
    arm_ids = sort!(unique(arms))
    if any(arm_ids .!= 0:length(arm_ids)-1)
        throw("Write code to fix this.")
    end
    study_ids = sort!(unique(studies))
    if any(study_ids .!= 0:length(studies_ids)-1)
        throw("Write code to fix this.")
    end
    decr = 0
    for s in 1:last(studies)
        
    end
end

abstract type AbstractNetworkEffect{P, T, L} <: PaddedMatrices.AbstractMutableFixedSizeVector{P, T, L} end
@inline Base.pointer(ne::AbstractNetworkEffect) = pointer(ne.effect)

abstract type AbstractFixedEffect{S, P, T, L} <: AbstractNetworkEffect{P, T, L} end
abstract type AbstractRandomEffect{P, T, L} <: AbstractNetworkEffect{P, T, L} end
struct FixedEffectM{S, P, T, L} <: AbstractFixedEffect{S, P, T, L}
    effect::FixedSizeVector{P,T,L}
    μ::FixedSizeVector{S,T,S}
    σᵦ::T
end
struct FixedEffect{S, P, T, L} <: AbstractFixedEffect{S, P, T, L}
    effect::PtrVector{P,T,L,false}
    μ::PtrVector{S,T,S}
    σᵦ::T
end
struct RandomEffectM{P, T, L} <: AbstractRandomEffect{P, T, L}
    effect::FixedSizeVector{P,T,L}
    σₑ::T
    σᵦ::T
end
struct RandomEffect{P, T, L} <: AbstractRandomEffect{P, T, L}
    effect::PtrVector{P,T,L,false}
    σₑ::T
    σᵦ::T
end

function ragged_network_meta_analysis_quote(
    israndom::Vector{Bool}, nmodelparams::Int, ntreatments::Int, (efitp, diitp, tritp)::NTuple{3,Bool},
    transform_v::Vector{Symbol}, arities::Vector{Int}, T, sptr::Bool, partial::Bool
)
    q = quote target = zero($T) end
    if partial
        diff_partials
    end
end

function gather_network_meta_analysis_quote(
    israndom::Vector{Bool}, nmodelparams::Int, ntreatments::Int, (efitp, diitp, tritp)::NTuple{3,Bool},
    transform_v::Vector{Symbol}, arities::Vector{Int}, T, sptr::Bool, partial::Bool 
)
    q = quote target = zero($T) end
    
end

const SUPPORTED_TRANSFORMS = Dict{Symbol,Int}(
    :ITP => 2,
    :Eₘₐₓ => 2,
    :exp => 1,
    :identity => 1,
    :* => 1,
    :logistic => 1
)
func_type_to_symbol(@nospecialize(f::Type{<:Function})) =  Symbol(string(f)[8:end-1])

"""

"""
function network_meta_analysis_quote(effects, differences, network, transforms, sptr::Bool = true, partial::Bool = true)
    effects_is_tuple = effects <: Tuple
    if effects_is_tuple
        israndom = Bool[ e <: AbstractRandomEffect for e in effects.parameters ]
    else
        israndom = Bool[ e <: AbstractRandomEffect ]
    end
    neffects = length(israndom)
    differences_is_tuple = differences <: Tuple
    nmodelparams = differences_is_tuple ? length(differences.parameters) : 1
    # Get the first two parameters of the MutableFixedSizeArray, they are Tuple{dims...}, eltype
    ntreatments_tup, T = (differences_is_tuple ? first(differences.parameters) : differences).parameters
    # Assume that length(ntreatments_tup.parameters) == 1, and extract the first element
    ntreatments = first(ntreatments_tup.parameters)::Int
    isragged = network <: RaggedNetwork
    transforms_is_tuple = false
    transform_v = if transforms === nothing
        [:identity]
    elseif transforms <: Tuple
        transforms_is_tuple = true
        [func_type_to_symbol(t) for t in transforms.parameters]
    else
        [func_type_to_symbol(transforms)]
    end
    @assert length(transform_v) == neffects
    arities = getindex.(SUPPORTED_TRANSFORMS, transform_v)
    @assert arities == nmodelparams
    tuple_checks = (effects_is_tuple, differences_is_tuple, transforms_is_tuple)
    if isragged
        ragged_network_meta_analysis_quote(israndom, nmodelparams, ntreatments, tuple_checks, transform_v, arities, T, sptr, partial)
    else
        gather_network_meta_analysis_quote(israndom, nmodelparams, ntreatments, tuple_checks, transform_v, arities, T, sptr, partial)
    end
end


"""
α = FixedEffect(  baselines, δα )
θ = RandomEffect( effects, δθ, stdev )
(α,θ) ~ NetworkMetaAnalysis( network, transforms )

How to pass in additional args?
"""
@generated function NetworkMetaAnalysis(sptr::StackPointer, effects::E, δ::D, network::N) where {E,D,N}
    network_meta_analysis_quote(E, D, N, nothing, sptr = true)
end
@generated function NetworkMetaAnalysis(sptr::StackPointer, effects::E, δ::D, intransforms::IT, network::N) where {E,D,N}
    network_meta_analysis_quote(E, D, N, nothing, sptr = true)
end
@generated function NetworkMetaAnalysis(sptr::StackPointer, effects::E, δ::D, network::N, outtransforms::OT) where {E,D,N,OT}
    network_meta_analysis_quote(E, D, N, nothing, sptr = true)
end
@generated function NetworkMetaAnalysis(sptr::StackPointer, effects::E, δ::D, intransforms::IT, network::N, outtransforms::OT) where {E,D,N,IT,OT}
    network_meta_analysis_quote(E, D, N, T, sptr = true)
end

@def_stackpointer_fallback NetworkMetaAnalysis ∂NetworkMetaAnalysis
function __init__()
    @add_stackpointer_method NetworkMetaAnalysis ∂NetworkMetaAnalysis
end

end # module
