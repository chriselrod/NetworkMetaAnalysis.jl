

"""
Calculates effects from differences.
"""    
function gather_fixed_effects(R, ptrsym, baseptr, basesym, diffsym, effectsym, partial = false, maskfirst = false)
    diff_ = Symbol(diffsym, :_)
    effect_ = Symbol(effectsym, :_)
    q = if maskfirst
        quote
            $(Symbol(effect_,1)) = vload(Vec{$W64,Float64}, $baseptr, $(VectorizationBase.max_mask(Float64)) ⊻ (one(MASK_TYPE) << go.offset - one(MASK_TYPE)))
        end
    else
        quote
            $(Symbol(effect_,1)) = vload(Vec{$W64,Float64}, $baseptr)
        end
    end
    q2 = quote
        $basesym = vsub( $(Symbol(effect_,1)), $(Symbol(diff_,1)) )
        baseptr += 8*$W64
    end
    append!(q.args, q2.args)
    for r ∈ 2:R
        push!(q.args, :($(Symbol(effect_, r)) = vadd($(Symbol(diff_, r)), $basesym)))
    end
    q
end


