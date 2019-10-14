abstract type AbstractFixedEffect{S, P, T, L} <: AbstractNetworkEffect{P, T, L} end
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




