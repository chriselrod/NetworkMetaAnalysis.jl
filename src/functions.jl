@inline function normal_lpdf(y, μ, τ, logτ)
    z = vsub(y, μ)
    vfmadd(τ, vmul(z,z), logτ)
end

@inline function ∂normal_lpdf(y, μ, τ, logτ, σ²)
    z = vsub(y, μ)
    zτ = vmul(z,τ)
    t = vfmadd(z, zτ, logτ)
    dtdy = zτ
    dtdμ = vsub(zτ)
    dtdτ = vfmadd(z, z, σ²)
    t, dtdy, dtdμ, dtdτ
end
