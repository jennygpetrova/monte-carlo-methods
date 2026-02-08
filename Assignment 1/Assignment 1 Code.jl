using Random, Primes
using Pkg
using StatsPlots
using Plots
using LaTeXStrings
using GLMakie
using Plots
using GLMakie
GLMakie.activate!()

# ---------------------
# Chapter 1 Problem 6
# ---------------------

# ------- 6 (a) -------

# test function for RANDU with seed x0=1
function randutest(n,x=1)
    M = 2^(31)
    for i in 1:n
        x = mod(65539 * x, M) # RANDU generator
    end
    return x
end
x = randutest(100)
println(x)

# complete RANDU generator
function randu(n,x=1)
    seq = zeros(n) 
    u_lower = [] 
    u_upper = [] 
    M = 2^(31)
    for i in 1:n
        seq[i] = x # store generated points
        x = mod(65539 * x, M) # RANDU generator
    end
    u = seq / M
    for i in 1:n
        # if 0.5 <= u_{i+1} <= 0.51, then store u_{i} and u_{i+2} as coordinate pairs
        if u[i] <= 0.51 && u[i] >= 0.50 
            push!(u_lower, u[i-1])
            push!(u_upper, u[i])
        end
    end
    return u, u_lower, u_upper
end

u1, u_lower, u_upper = randu(20002)

# create 2D scatterplot of u_{i} vs u_{i+2}
Plots.scatter(u_lower, u_upper,
    xlabel = L"$u_{i}$",
    ylabel = L"$u_{i+2}$",
    title = "RANDU",
    legend = false
)


# ------- 6 (b) -------

u2, _, _ = randu(1002)

# store consecutive triplets
x = u2[1:1000]
y = u2[2:1001]
z = u2[3:1002]

# create 3D (interactive) scatterplot of triplets (u_{i} vs u_{i+1} vs u_{i+2})
fig = Figure()
ax = Axis3(fig[1, 1], xlabel = L"$u_i$", ylabel = L"$u_{i+1}$", zlabel = L"$u_{i+2}$")
GLMakie.scatter!(ax, x, y, z; markersize = 2)
display(fig)


# ---------------------
# Chapter 1 Problem 10
# ---------------------

m = 40
N1 = 1000
N2 = 10000
s = 10

# create generator for MC method
rng = MersenneTwister(1234)

# create generator for RQMC method
lprimes=primes(2,200)

# generate nth element of the van der Corput sequence in base b
function vdc(n,b)
    digs=digits(n,base=b)
    sum=0.0
    for i in 1:length(digs)
        sum = sum + digs[i]/b^i
    end
    return(sum)
end

# get the nth element of the s-dimensional Halton sequence
function halton(n,s)
    [vdc(n,i) for i in lprimes[1:s]]
end

# randomly shift the s-dimensional Halton sequence
function shifthalton(n,s)
    u = rand(rng, s)
    map(k->mod(k,1),halton(n,s).+u)
end

# estimate the integral for Nxs matrix with N pseudorandomly generated s-dimensional sequences
function fintegral(N, s, method="rqmc")
    I = 0.0
    for i in 1:N
        if method=="mc"
            q = rand(rng, s)
        else
            q = shifthalton(i, s)
        end
        term = 1.0
        for j in 1:(length(q))
            term *= (pi/2) * sin(pi * q[j])
        end
        I += term
    end
    return (I/N)
end

# print integral estimates for MC and RQMC, for N=1000 and N=10,000
println("Integral Estimates for N=$N1")
I_mc1 = [fintegral(N1, s, "mc") for _ in 1:m]
println("Mersenne Twister: $I_mc1")
I_rqmc1 = [fintegral(N1, s) for _ in 1:m]
println("Randomly shifted Halton: $I_rqmc1")
println("Integral Estimates for N=$N2")
I_mc2 = [fintegral(N2, s, "mc") for _ in 1:m]
println("Mersenne Twister: $I_mc2")
I_rqmc2 = [fintegral(N2, s) for _ in 1:m]
println("Randomly shifted Halton: $I_rqmc2")

# create box-and-whisker plots for each set of integral estimates
data=[I_mc1, I_rqmc1, I_mc2, I_rqmc2]
Plots.boxplot(data, legend=false)
Plots.xticks!(1:4, ["MC (N=1K)", "RQMC (N=1K)", "MC (N=10K)", "RQMC (N=10K)"])
Plots.ylims!(0.7,1.3)