### Integration methods ###
#const figPATH = "/Users/pau/Desktop/MSc_CNS/Thesis/Figures"
#const netsimPATH = "/Users/pau/Desktop/MSc_CNS/Thesis/Code/networks/simulations"
using Pkg

function check_pkg(pkg::AbstractString)
    print("install uninstalled packages? (y/n) ")
    n = readline()
    if n == "y"
        try
            pkg != "METADATA" && pkg != "REQUIRE" && pkg[1] != '.' && Pkg.cd(isdir, pkg)
        catch
            Pkg.add(pkg)
        end
    else
        println("packages are needed to continue.")
        return nothing
    end
end
function check_pkg(pkgs::Vector{String})
    println("install uninstalled packages? (y/n)")
    n = readline()
    if n == "y"
        for pkg in pkgs
            try
                pkg != "METADATA" && pkg != "REQUIRE" && pkg[1] != '.' && Pkg.cd(isdir, pkg)
            catch
                Pkg.add(pkg)
            end
        end
    else
        println("packages are needed to continue.")
        return nothing
    end
end

function Euler( g::Function, vars::Vector{Float64}, params::Vector{Float64}
              , h::Float64)
    vars .+= g(vars, params).*h
    return vars
end

function Euler( g::Function, var::Float64, params::Vector{Float64}
              , h::Float64)
    return g(var, params)*h
end

function Heuns( f::Function, vars::Vector{Float64}, params::Vector{Float64}
              , h::Float64)
    dvars = f(vars, params)
    int_vars = vars .+ h .* dvars
    int_dvars = f(int_vars, params)
    @. vars += h * 0.5 * (dvars + int_dvars)
    return vars
end

function Heuns2( f::Function, vars::Vector{Float64}, params::Vector{Float64}
              , h::Float64)
    return h * 0.5 .* (dvars .+ f(vars .+ h .* f(vars, params), params))
end

"""
    4th order Runge Kuta method for a system of ODEs
"""
function RK4( f::Function, vars::Vector{Float64}, params::Vector{Float64}
            , h::Float64)

    k1 = f(vars, params)
    vars_k2 = vars .+ 0.5 * h .* k1
    k2 = f(vars_k2, params)
    vars_k3 = vars .+ 0.5 * h .* k2
    k3 = f(vars_k3, params)
    vars_k4 = vars .+ h .* k3
    k4 = f(vars_k4, params)

    slope = (k1 .+ 2.0.*k2 .+ 2.0.*k3 .+ k4) ./ 6.0
    @. vars = vars + (h * slope)
    return vars
end

### END Integration methods ###

function overprint(str; i=1)
    if i==1
        println(str)
    elseif i>1
        print("\u1b[1F")
        #Moves cursor to beginning of the line n (default 1) lines up
        print(str)   #prints the new line
        print("\u1b[0K")


        # clears  part of the line.
        #If n is 0 (or missing), clear from cursor to the end of the line.
        #If n is 1, clear from cursor to beginning of the line.
        #If n is 2, clear entire line.
        #Cursor position does not change.

        println() #prints a new line, i really don't like this arcane codes
    end

end
