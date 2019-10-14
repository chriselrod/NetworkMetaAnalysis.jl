

"""
Quote accumulates separate target and target_base.

Calculates differences from effects, to apply normal distribution to them.
"""
function gather_random_effects(R, ptrsym, effectptr, effectsym, diffsym, partial = false, maskfirst = false)
    effect_ = Symbol(effectsym, :_)
    diff_ = Symbol(diffsym, :_)
    q = if maskfirst
        quote
            $(Symbol(effect_,1)) = vload(Vec{$W64,Float64}, $effectptr, $(VectorizationBase.max_mask(Float64)) ⊻ (one(MASK_TYPE) << go.offset - one(MASK_TYPE)))
            effectptr += 8*$W64
            $([:(($(Symbol(effect_,r)) = vload(Vec{$W64,Float64},$effectptr)); $effectptr += 8*$W64) for r ∈ 2:R]...)
        end
    else
        quote
            $([:(($(Symbol(effect_,r)) = vload(Vec{$W64,Float64},$effectptr)); $effectptr += 8*$W64) for r ∈ 1:R]...)
        end
    end
    effect_1 = Symbol(effect_,1)
    q2 = quote
        b = vsub($effect_1, $(Symbol(diff_,1)))
        s = vbroadcast(Vec{$W64,Float64}, 0.0)
        target_base = vfmadd($effect_1, $effect_1, target_base)
    end
    append!(q.args, q2.args)
    for r ∈ 2:R
        qr = quote
            diff = vsub(vsub($(Symbol(effect_,r)),$(Symbol(diff_,r))), b)
            target = vadd(
                normal_lpdf(
                    diff, $(r == 2 ? :s : :(vmul(vbroadcast(Vec{$W64,Float64}, $(1/(r-1))), s))), $(Symbol(:τ_,r)), $(Symbol(:logrootτ_,r))
                ), target
            )
        end
        r < R && push!(qr.args, :(s = vadd(s, diff)))
        append!(q.args, qr.args)
    end
    q
end


