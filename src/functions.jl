
const VectorOrDouble{T} = Union{T,<:PaddedMatrices.AbstractFixedSizeVector{<:Any,T}}

struct Eₘₐₓ{T,N<:AbstractNetwork,M<:VectorOrDouble{T},D<:VectorOrDouble{T},B<:Union{T,Nothing}}
    N::N
    Eₘₐₓ::M
    ED₅₀::D
    μ::B
end

@inline Eₘₐₓ(em::M, ed::D) where {M,D} = Eₘₐₓ{M,D,Nothing}(em, ed, nothing)

@generated function ProbabilityDistributions.Normal(
    data::AbstractVector{T},
    emax::Eₘₐₓ{T,N,M,D,B},
    τ::AbstractVector{T},
    ::Val{track} = Val((false,true,false))
) where {T,R,N<:GatherNetwork{R},M,D,B,track}
    track_data, track_emax, track_τ = track
    @assert !(track_data | track_τ) # for now; ought to implement later.

    q = quote
        $(Expr(:meta,:inline)) # we need to inline to avoid allocations from actually constructing Eₘₐₓ object
        
    end
    
    
end

@generated function ProbabilityDistributions.∂Normal(
    sp::StackPointer,
    data::AbstractVector{T},
    emax::Eₘₐₓ{T,N,M,D,B},
    τ::AbstractVector{T},
    ::Val{track} = Val((false,true,false))
) where {T,R,N<:GatherNetwork{R},M,D,B,track}

    track_data, track_emax, track_τ = track
    @assert !(track_data | track_τ) # for now; ought to implement later.

    q = quote
        $(Expr(:meta,:inline)) # we need to inline to avoid allocations from actually constructing Eₘₐₓ object
        
    end
    

end



