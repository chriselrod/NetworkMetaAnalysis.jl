
struct Eₘₐₓ

end

@genereated function ProbabilityDistributions.Normal(
    data::AbstractVector{T},
    μ::T
    emax::Eₘₐₓ,
    τ::AbstractVector{T},
    ::Val{track} = Val((false,true,false))
) where {T,track}
    track_data, track_emax, track_τ = track
    @assert !(track_data | track_τ) # for now; ought to implement later.

    
    
end


