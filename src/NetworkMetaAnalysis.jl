module NetworkMetaAnalysis

using VectorizationBase, SIMDPirates, LoopVectorization,
    StackPointers, PaddedMatrices, ReverseDiffExpressionsBase

export NetworkMetaAnalysis, RaggedNetwork, GatherNetwork,
    FixedEffectParameters, FixedEffect,
    MetaAnalysis, ∂MetaAnalysis

abstract type AbstractNetwork end
abstract type AbstractNetworkEffect{T} <: DenseVector{T} end
# abstract type AbstractNetworkEffect{P, T, L} <: PaddedMatrices.AbstractMutableFixedSizeVector{P, T, L} end
# @inline Base.pointer(ne::AbstractNetworkEffect) = pointer(ne.effect)

const W64, Wshift64 = VectorizationBase.pick_vector_width_shift(Float64)
const MASK_TYPE = VectorizationBase.mask_type(Float64) # UInt8, unless we get something like avx1024.

# @inline smask(r) = VectorizationBase.max_mask(Float64) ⊻ ( (one(MASK_TYPE) << (r & MASK_TYPE(8sizeof(MASK_TYPE)-1))) - one(MASK_TYPE))

"""
Represents a network as a ragged matrix.
"""
struct RaggedNetwork <: AbstractNetwork
    
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
function network_meta_analysis_quote(effects, differences, network, transforms; sptr::Bool = true, partial::Bool = false)
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
    # ntreatments_tup, T = (differences_is_tuple ? first(differences.parameters) : differences).parameters
    ntreatments_tup = first((differences_is_tuple ? first(differences.parameters) : differences).parameters)
    # Assume that length(ntreatments_tup.parameters) == 1, and extract the first element
    ntreatments = first(ntreatments_tup.parameters)::Int
    isragged = network <: RaggedNetwork
    transforms_is_tuple = false
    transform_v = if transforms === nothing
        [:identity for _ ∈ 1:differences]
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
        ragged_network_meta_analysis_quote(israndom, nmodelparams, ntreatments, tuple_checks, transform_v, arities, sptr, partial)
    else
        arm_lengths = [Int(r) for r ∈ first(network.parameters)]
        gather_network_meta_analysis_quote(israndom, nmodelparams, ntreatments, arm_lengths, tuple_checks, transform_v, arities, sptr, partial)
    end
end


"""
α = FixedEffect(  baselines, δα )
θ = RandomEffect( effects, δθ, stdev )
(α,θ) ~ NetworkMetaAnalysis( network, transforms )

How to pass in additional args?
"""
@generated function MetaAnalysis(sptr::StackPointer, effects::E, δ::D, network::N) where {E,D,N}
    network_meta_analysis_quote(E, D, N, nothing, sptr = true)
end
@generated function MetaAnalysis(sptr::StackPointer, effects::E, δ::D, intransforms::IT, network::N) where {E,D,IT,N}
    network_meta_analysis_quote(E, D, N, nothing, sptr = true)
end
@generated function MetaAnalysis(sptr::StackPointer, effects::E, δ::D, network::N, outtransforms::OT) where {E,D,N <: AbstractNetwork,OT}
    network_meta_analysis_quote(E, D, N, nothing, sptr = true)
end
@generated function MetaAnalysis(sptr::StackPointer, effects::E, δ::D, intransforms::IT, network::N, outtransforms::OT) where {E,D,N,IT,OT}
    network_meta_analysis_quote(E, D, N, T, sptr = true)
end

function ∂MetaAnalysis end

include("functions.jl")
include("fixed_effects.jl")
include("gather_fixed_effects.jl")
# include("random_effects.jl")
# include("gather_random_effects.jl")
include("gather_network.jl")

# include("ragged_fixed_effects.jl")
# include("ragged_random_effects.jl")
# include("ragged_network.jl")


@def_stackpointer_fallback MetaAnalysis ∂MetaAnalysis
function __init__()
    @add_stackpointer_method MetaAnalysis ∂MetaAnalysis
end

end # module
