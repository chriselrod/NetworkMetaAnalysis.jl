

"""
Calculates effects from differences.
"""    
function gather_fixed_effects(R, ptrsym, baseptr, basesym, diffsym, effectsym, partial = false, masklast = false, S)
    diff_ = Symbol(diffsym, :_)
    effect_ = Symbol(effectsym, :_)
    q = if masklast
        quote
            $(Symbol(effect_,1)) = vload(Vec{$W64,Float64}, $baseptr, $(VectorizationBase.mask(Float64, S & (W64-1))))
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
        push!(q.args, :( $(Symbol(effect_, r)) = vadd($(Symbol(diff_, r)), $basesym) ))
    end
    q
end

function store_fixed_effect(R, ptrfixed, effectsym, maskfirst = false)
    effect_ = Symbol(effectsym, :_)
    q = if maskfirst
        quote
            vstore!( $ptrfixed, $(Symbol(effect_,1)), maskf)
            $ptrfixed += 8 * maskfcount
        end
    else
        quote
            vstore!( $ptrfixed, $(Symbol(effect_,1)) ); $ptrfixed += 8*$W64
        end
    end
    push!(q.args, :(vstore!($ptrfixed, $(Symbol(effect_,r)), $(Symbol(:mask_,r))); $ptrfixed += 8*$W64))
    for r ∈ 2:R
        push!(q.args, :(vstore!($ptrfixed, $(Symbol(effect_,r)), $(Symbol(:mask_,r))); $ptrfixed += $(Symbol(:maskcount_,r))))
    end
    q
end



