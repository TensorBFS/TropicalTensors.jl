using .CuYao

function _init_reg(::Type{T}, L::Int, usecuda::Val{:true}) where T
    reg = ArrayReg(CuArrays.ones(Tropical{T}, 1<<L))
end
