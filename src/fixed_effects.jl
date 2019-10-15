
# abstract type AbstractFixedEffect{S, P, T, L} <: AbstractNetworkEffect{P, T, L} end

struct FixedEffect{T, V <: DenseVector{T}} <: AbstractNetworkEffect{T} # AbstractFixedEffect{P, T, L}
    data::V#ector{T}
end
Base.length(fe::FixedEffect) = length(fe.data)
Base.size(fe::FixedEffect) = size(fe.data)
Base.size(fe::FixedEffect, n) = size(fe.data, n)
Base.strides(fe::FixedEffect) = strides(fe.data)
Base.stride(fe::FixedEffect, n) = stride(fe.data, n)
Base.pointer(fe::FixedEffect) = pointer(fe.data)
Base.getindex(fe::FixedEffect, I...) = getindex(fe.data, I...)
Base.setindex!(fe::FixedEffect, v, I...) = setindex!(fe.data, v, I...)
Base.IndexStyle(::Type{<:FixedEffect}) = IndexLinear()

struct FixedEffectParameters{D,S,T}
    μ::PtrVector{S,T,S,false}
    δ::PtrVector{D,T,D,false}
end





