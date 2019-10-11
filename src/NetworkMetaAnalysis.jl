module NetworkMetaAnalysis

using VectorizationBase, SIMDPirates,
    StackPointers, PaddedMatrices, ReverseDiffExpressionsBase

export NetworkMetaAnalysis, RaggedNetwork, GatherNetwork

abstract type AbstractNetwork end
abstract type AbstractNetworkEffect{P, T, L} <: PaddedMatrices.AbstractMutableFixedSizeVector{P, T, L} end
@inline Base.pointer(ne::AbstractNetworkEffect) = pointer(ne.effect)

const W64, Wshift64 = VectorizationBase.pick_vector_width_shift(Float64)
const MASK_TYPE = VectorizationBase.mask_type(Float64) # UInt8, unless we get something like avx1024.

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

Base.Sort.defalg(::SortingSet) = Base.Sort.QuickSortAlg()


struct GatherOffsets
    offsets::Matrix{Int}#Vector{NTuple{W64,Core.VecElement{Int}}}
    masks::Vector{MASK_TYPE}
    nreps::Vector{Int}
    offset::Int # = rem - W
    # increment::Vector{Int}
end

function cumsum0(x::AbstractArray{T}) where {T}
    y = similar(x)
    y[1] = zero(T)
    s = zero(T)
    for i ∈ 2:length(x)
        y[i] = s += x[i-1]
    end
    y
end

# """
# GatherOffsets( items, counts )
 # - Items, the offsets of the things we wish to load into registers. For example, offsets of arm differences.
 # - Lanes, the category we wish to divide things by. For example, studies, where we wish to load each arm corresponding to that study in the same lane.
 # - Counts, the number of items corresponding to each lane. For example, the number of arms per study.

# counts are counts of the things you wish to load into registers.
# """
function GatherOffsets(items, counts, csc = cumsum0(counts))
    # @show items'
    # @show counts'
    # @show csc
    ncounts = length(counts)
    count_rem = ncounts & (W64-1)
    count_reps = ncounts >>> Wshift64

    count_has_remainder = count_rem != 0
    if count_has_remainder
        count_per_lane = Vector{Int}(undef, count_reps + 1)
        dim2 = counts[count_rem]
        count_per_lane[1] = dim2
        # offset = count_rem - W64
    else
        count_per_lane = Vector{Int}(undef, count_reps)
        dim2 = 0
        # offset = 0
    end
    offset = count_rem - W64
    for i ∈ 1:count_reps
        c = counts[count_rem + (i << Wshift64)]
        dim2 += c
        count_per_lane[i + count_has_remainder] = c
    end
    
    Vf64 = Vec{W64,Float64}
    Vint = Vec{W64,Int}

    # increment = Vector{Int}(undef, dim2+1); increment[1] = offset
    # offsets = Vector{NTuple{W64,Core.VecElement{Int}}}(undef, dim2) # gives offsets for vectorizing studies
    offsets = fill(-Int(999), W64, dim2) # May not need mask vector? mask = offsets >= 0
    masks = Vector{MASK_TYPE}(undef, dim2)
    ptr_counts = pointer(counts)
    sizeof_I = sizeof(eltype(counts))
    if count_has_remainder
        ptr_counts += offset*sizeof_I
    end
    offset_ptr = reinterpret(Ptr{Int}, pointer(offsets))
    ind = 1
    # @show !count_has_remainder, count_reps
    # @show offset, length(items)
    for s ∈ (!count_has_remainder):count_reps
        if s == 0
            count_inds = ntuple(w -> (w + offset) > 0 ? csc[w + offset] : -99999, Val(W64))
            mask = VectorizationBase.max_mask(Float64) ⊻ ((one(MASK_TYPE) << (W64-1)) - one(MASK_TYPE))
            count_vec = vload( Vint, ptr_counts, mask )
        else
            count_inds = ntuple(w -> csc[w + offset + (s << Wshift64)], Val(W64))
            count_vec = vload( Vint, ptr_counts )
        end
        # @show count_inds
        maxcount = count_per_lane[s + count_has_remainder]
        # @show count_vec
        for a ∈ 1:maxcount
            # @show a, ind, maxcount, s
            mask = SIMDPirates.vgreater_or_equal_mask( count_vec, vbroadcast(Vint, a) )
            # @show bitstring(mask)
            # Fill in offset_α
            for w ∈ 0:W64-1
                flag = one(MASK_TYPE) << w
                if (mask & flag) != zero(MASK_TYPE) # arm present
                    i = count_inds[w+1]
                    # @show i
                    count_inds = Base.setindex(count_inds, i+1, w+1)
                    itemsᵢ = items[i+1]
                    # Sentinel value; this is a baseline / something we do not load -> set corresponding mask to 0
                    if itemsᵢ < 0
                        mask ⊻= flag
                    else
                        VectorizationBase.store!( offset_ptr + w*sizeof(Int), itemsᵢ )
                    end
                end
            end
            masks[ind] = mask
            offset_ptr += W64*sizeof(Int)
            ind += 1
            # increment[ind] = count_ones(mask) # instruction fast enough that I should use it instead of loading.
        end
        ptr_counts += W64 * sizeof_I
    end
    GatherOffsets( offsets, masks, count_per_lane, offset ) # increment )
end

 
"""
Represents a network as a series of offsets.
"""
struct GatherNetwork <: AbstractNetwork
    α::GatherOffsets
    δ::GatherOffsets
end

function count_unique(v::AbstractVector{T}) where {T}
    d = Dict{T,Int}()
    for vᵢ ∈ v
        d[vᵢ] = get(d, vᵢ, 0) + 1
    end
    d
end
function count_order_mapping(d::Dict{K,V}, baseline = 1; rev = false) where {K,V}
    sign = rev ? -1 : 1
    vk = sort!([(v,k) for (k,v) ∈ d], lt = (s,a) -> isless((sign*s[1],s[2]),(sign*a[1],a[2])))
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
    arms_per_study = count_unique(studies)
    arm_frequency = count_unique(arms)
    ss = sort!(
        SortingSet(copy(studies), copy(arms), copy.(doses)),
        lt = (s,a) -> isless( # More frequently included treatments (arm freq) go first 
            (arms_per_study[s[1]],s[1],-arm_frequency[s[2]],s[2]),
            (arms_per_study[a[1]],a[1],-arm_frequency[a[2]],a[2])
        )
    )
    # @show arm_frequency
    study_map, study_count = count_order_mapping(arms_per_study, 0, rev = false)
    arm_map, arm_count = count_order_mapping(arm_frequency, -1, rev = true)

    # @show typeof(ss)
    # @show typeof(ss.x), typeof(study_map)
    studies = [study_map[x] for x ∈ ss.x]
    arms = [arm_map[y] for y ∈ ss.y]
    # studies and arms now have the correct ids / names.
    # @show studies'
    # @show arms'
    # @show count_unique(arms)
    study_offsets = GatherOffsets( arms, study_count )
    # Need to construct study_offset_vec like the arms above; -1 indicates skipping.
    
    ad = [ sizehint!(Int[], c) for c ∈ arm_count ]
    # @show study_offsets.offsets'
    # @show arm_count'
    for k ∈ 1:length(study_offsets.offsets)
        o = study_offsets.offsets[k]
        o < 0 && continue
        # o gives the arm
        push!(ad[o+1], k-1)
    end
    adp = vcat(ad...)
    # @show extrema(adp)
    arm_offsets = GatherOffsets( adp, arm_count[2:end] )
    # arm_offsets = GatherOffsets( vcat(ad...), arm_count )

    GatherNetwork( study_offsets, arm_offsets ), ss.z
end


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

@def_stackpointer_fallback MetaAnalysis ∂MetaAnalysis
function __init__()
    @add_stackpointer_method MetaAnalysis ∂MetaAnalysis
end

end # module
