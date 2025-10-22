using Statistics
using Plots
pyplot()
function test(n)
    x = zeros(n)
    for i=1:n
        x[i] = round(Int, N_start + Nstartsd*randn())
        #x[i]   = round(Int, N_start + Nstartsd*randn() + 0.5)
    end
    println()
    println(mean(x))
    println(std(x))
    return x
end
x = test(10000000)
histogram(x)
