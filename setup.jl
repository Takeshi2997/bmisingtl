module Const

    struct Param

        # System Size
        l::Int64 # lattice length
        n::Int64 # site number
        x::Array{Int64, 2} # coordinate

        # Learning Constant
        iters_num::Int64 # Gibbs step number
        burnintime::Int64 # burn in time
        sample_interval::Int64 # sample interval
    end

   l = 8
   n = l^2
   x = reshape(1:n, l, l)

   iters_num = 50000
   burnintime = 1000
   sample_interval = 1000
end
