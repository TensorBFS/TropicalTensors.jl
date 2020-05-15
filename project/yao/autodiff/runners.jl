using TropicalTensors

function benchmarker(L)
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    A = stack4reg(reg, L)
    B = stack4reg(reg, L-1)
    isolve(Float32(0.0), reg, L, ones(Float32, L*(L-1)*2), A, B)
    println("cached array size: $(count(x->x isa AbstractArray, NiLang.GLOBAL_STACK))")
    empty!(NiLang.GLOBAL_STACK)
end

function benchmarker1(L)
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    isolve_largemem(Float32(0.0), reg, L, ones(Float32, L*(L-1)*2))
    println("cached array size: $(count(x->x isa AbstractArray, NiLang.GLOBAL_STACK))")
    empty!(NiLang.GLOBAL_STACK)
end

function benchmarker2(L)
    reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
    _spinglass_yao(reg, L, ones(Float32, L*(L-1)*2))
end

function benchmarker4(L)
    reg = ArrayReg(Tropical.(map(x->ForwardDiff.Dual(x, zero(x)), zeros(1<<L))))
    Js = ones(Float32, L*(L-1)*2)
    dJs = map(i->ForwardDiff.Dual(Js[i], i==1 ? one(Js[i]) : zero(Js[i])), 1:length(Js))
    _spinglass_yao(reg, L, dJs)
end

# fail
function benchmarker5(L)
    function f(Js)
        reg = ArrayReg(ones(Tropical{Float32}, 1<<L))
        _spinglass_yao(reg, L, Js).n
    end
    Zygote.gradient(f, ones(Float32, L*(L-1)*2))
end
