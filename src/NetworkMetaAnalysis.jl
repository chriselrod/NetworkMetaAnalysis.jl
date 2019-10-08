module NetworkMetaAnalysis

using VectorizationBase, SIMDPirates,
    StackPointers, PaddedMatrices, ReverseDiffExpressionsBase

export NetworkMetaAnalysis, RaggedNetwork, GatherNetwork

abstract type AbstractNetwork end

"""
Represents a network as a ragged matrix.
"""
struct RaggedNetwork <: AbstractNetwork

end

"""
Used for sorting sets of values together without using sortperm.
"""
struct SortingSet{N,T} <: AbstractVector{Tuple{Int,Int,NTuple{N,T}}}
    x::Vector{Int}
    y::Vector{Int}
    z::NTuple{N,Vector{T}}
end
Base.size(s::SortingSet) = size(s.x)
Base.length(s::SortingSet) = length(s.x)
Base.IndexStyle(::Type{<:SortingSet}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(s::SortingSet{N}, i) where {N}
    (s.x[i], s.y[i], ntuple(n -> s.z[n][i], Val(N)))
end
Base.@propagate_inbounds function Base.setindex!(
    s::SortingSet{N,T}, (x,y,z)::Tuple{Int,Int,NTuple{N,T}}, i
) where {N,T}
    (setindex!(s.x, x, i), setindex!(s.y, y, i), ntuple(n -> setindex!(s.z[n], z[n], i), Val(N)))
end

function sortpermsort(ss::SortingSet{N}) where {N}
    sp = sortperm(ss.x)
    SortingSet(
        ss.x[sp], ss.y[sp],
        ntuple(n -> ss.z[n][sp], Val(N))
    )
end
Base.Sort.defalg(::SortingSet) = Base.Sort.QuickSortAlg()


"""
Represents a network as a series of offsets.
"""
struct GatherNetwork{U} <: AbstractNetwork
    arms_and_studies::Matrix{Int}
    masks::Vector{U}
    narms::Vector{Int}
end

function count_unique(v::AbstractVector{T}) where {T}
    d = Dict{T,Int}()
    for vᵢ ∈ v
        d[vᵢ] = get(d, vᵢ, 0) + 1
    end
    d
end
function count_order_mapping(d::Dict{K,V}, baseline = 1) where {K,V}
    vk = sort!([(v,k) for (k,v) ∈ d])
    m = Dict{K,Int}()
    count = Vector{V}(undef, length(vk))
    ind = baseline
    for (v,k) ∈ vk
        m[k] = ind
        ind += 1
        count[ind - baseline] = v
    end
    m, count
end

function GatherNetwork(studies, arms, doses...)
    nstudyarms = length(studies)
    @boundscheck begin
        all(d -> length(d) == nstudyarms, doses) || throw("There should be one dose per study ID.")
        nstudyarms == length(arms) || throw("There should be one study ID per arm.")
    end
    arms_per_study = count_unqiue(studies)
    arm_frequency = count_unique(arms)
    ss = sort(
        SortingSet(studies, arms, doses),
        lt = (s,a) -> isless( # More frequently included treatments (arm freq) go first 
            (arms_per_study[s[1]],s[1],-arm_frequency[s[2]]),
            (arms_per_study[a[1]],a[1],-arm_frequency[a[2]])
        )
    )
    study_map, study_count = count_order_mapping(arms_per_study, 1)
    arm_map, arm_count = count_order_mapping(arm_frequency, -1)

    studies = getindex.(study_map, ss.x)
    arms = getindex.(arm_map, ss.y)
    # studies and arms now have the correct ids / names.

    W, Wshift = VectorizationBase.pick_vector_width(Float64)
    
    # offsets = Array{Int}( undef, W, ( nstudyarms + W - 1 ) >> Wshift, 2 )
    # rem = nstudyarms & (W - 1)
    nstudies = length(study_count)
    study_rem = nstudies & (W-1)
    # study_slack = study_rem == 0 ? 0 : W - studyrem
    study_reps = nstudies >>> Wshift

    study_has_remainder = study_rem != 0
    if study_has_remainder
        arms_per_set = Vector{Int}(undef, study_rem + 1)
        dim2 = study_count[study_rem]
        arms_per_set[1] = dim2
        so = W - study_rem
    else
        arms_per_set = Vector{Int}(undef, study_rem)
        dim2 = 0
        so = 0
    end
    for i ∈ 1:study_reps
        sc = study_count[study_rem + (i << Wshift)]
        dim2 += sc
        arms_per_set[i + study_has_remainder] = sc
    end
    
    base_mask = VectorizationBase.max_mask(Float64)
    Vf64 = VectorizationBase.pick_vector(Float64)
    Vi64 = VectorizationBase.pick_vector(Int64)

    offset_α = Array{Int}(undef, W, dim2, 2) # gives offsets for vectorizing studies
    offset_δ = Array{Int}(undef, W, dim2) # gives offsets for vectorizing δs, so we can calculate ∂s from stored ∂out∂δ
    masks = Vector{typeof(base_mask)}(undef, dim2)
    ptr_counts = pointer(study_count)
    sizeof_I = sizeof(eltype(study_count))
    if study_has_remainder
        ptr_counts -= so*sizeof_I
    end
    ind = 1
    csc = cumsum(study_count)
    for s in (!study_has_remainder):study_reps
        if s == 0
            study_inds = ntuple(w -> (w - so) > 0 ? csc[w - so] : 0, pick_vector_width_val(Float64))
            mask = base_mask ⊻ ((one(base_mask) << (W-1)) - one(base_mask))
            narms_vec = vload( Vi64, ptr_counts, mask )
        else
            study_inds = ntuple(w -> csc[w - so + (s << Wshift)], pick_vector_width_val(Float64))
            narms_vec = vload( Vi64, ptr_counts )
        end
        narms = arms_per_set[s + study_has_remainder]
        for a ∈ 1:narms
            mask = SIMDPirates.vecbool_to_unsigned(
                SIMDPirates.vgreater_or_equal( narms_vec, vbroadcast(Vi64, a) )
            )
            masks[ind] = mask
            # Fill in offset_α
            for w in 1:W
                if (mask & (one(mask) << (w-1))) != zero(mask) # arm present
                    si = study_inds[w]
                    study_inds = setindex(study_inds, si + 1, w)
                    offset_α[w,ind,1] = studies[si]
                    offset_α[w,ind,2] = arms[si]
                end
            end            
            ind += 1
        end
        ptr_counts += W * sizeof(T)
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
