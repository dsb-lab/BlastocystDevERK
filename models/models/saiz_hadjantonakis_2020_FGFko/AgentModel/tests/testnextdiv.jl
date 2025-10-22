using Statistics
using Plots
pyplot()

function test(n)
    nextdiv1 = zeros(n)
    nextdiv2 = zeros(n)
    nor = Normal()
    ct=0
    for i=1:n
        nextdiv1[i] = ct + tdiv + (rand()*2*sdiv - sdiv)
        nextdiv2[i] = (1 - sdiv + rand()*2.0*sdiv)*tdiv
    end
    println(mean(nextdiv1))
    println(std(nextdiv1))
    println(mean(nextdiv2))
    println(std(nextdiv2))
    return nextdiv1, nextdiv2
end

nd1, nd2 = test(100000)
histogram(nd1)
histogram(nd2)
