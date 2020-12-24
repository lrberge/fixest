#=
Benchmark:
- Julia version: 1.2.0
- Author: laurent.berge
- Date: 2019-11-13
- ~: Replication in Julia
=#

#=
using Pkg
Pkg.add("RData")
Pkg.add("DelimitedFiles")
Pkg.add("FixedEffectModels")
Pkg.add("DataFrames")
Pkg.add("CSV")

pwd()
cd("./Google Drive/R_packages/fixest/fixest/_BENCHMARK")
=#

println("Loading packages...")
using RData, DelimitedFiles, FixedEffectModels, DataFrames, CSV
println("done.")

println("Loading data...")
base_all = load("DATA/base_all_simulations.RData", convert=true);
base_all = base_all["base_all"];
println(".")
base_10M = CSV.read("DATA/base_10M.csv", DataFrame)
println("done.")

# warming up
println("warming up...")
base = base_all[1]
test = @elapsed reg(base, @formula(ln_y~X1 + fe(dum_1)))
println("done.")

# 1 FE

# We run all the models
println("One FE.")
timings = Vector{Float64}(undef, 50)
for i=1:40
    println("i = ", i)
    base = base_all[i]
    print(size(base))
    timings[i] = @elapsed reg(base, @formula(ln_y~X1 + fe(dum_1)))
end

# The 10M DB
for i=1:10
    println("i = ", i + 40)
    timings[i + 40] = @elapsed reg(base_10M, @formula(ln_y~X1 + fe(dum_1)))
end

# saving
open("DATA/julia_bench_1FE.txt", "w") do io
           writedlm(io, timings)
       end

# 2 FE

# We run all the models
println("Two FEs.")
timings = Vector{Float64}(undef, 50)
for i=1:40
    println("i = ", i)
    base = base_all[i]
    timings[i] = @elapsed reg(base, @formula(ln_y~X1 + fe(dum_1) + fe(dum_2)))
end

# The 10M DB
for i=1:10
    println("i = ", i + 40)
    timings[i + 40] = @elapsed reg(base_10M, @formula(ln_y~X1 + fe(dum_1) + fe(dum_2)))
end

# saving
open("DATA/julia_bench_2FE.txt", "w") do io
           writedlm(io, timings)
       end


# 3 FE

# We run all the models
println("Three FEs.")
timings = Vector{Float64}(undef, 50)
for i=1:40
    println("i = ", i)
    base = base_all[i]
    timings[i] = @elapsed reg(base, @formula(ln_y~X1 + fe(dum_1) + fe(dum_2) + fe(dum_3)))
end

# The 10M DB
for i=1:10
    println("i = ", i + 40)
    timings[i + 40] = @elapsed reg(base_10M, @formula(ln_y~X1 + fe(dum_1) + fe(dum_2) + fe(dum_3)))
end

# saving
open("DATA/julia_bench_3FE.txt", "w") do io
           writedlm(io, timings)
       end



#
# The "difficult" data
#

println("Loading difficult data...")
base_diff = load("DATA/base_all_diff.RData", convert=true);
base_diff = base_diff["base_all_diff"];
println("done.")

# warming up
println("warming up...")
base = base_diff[1]
test = @elapsed reg(base, @formula(y ~ x1 + x2 + fe(id_indiv)))
println("done.")

timings = Vector{Float64}(undef, 120)

global index = 1
for size=1:4
    println("size = ", size)
    base = base_diff[size]
    for g=1:3
        println("g = ", g)
        for r=1:10
            if g == 1
                global timings[index] = @elapsed reg(base, @formula(y ~ x1 + x2 + fe(id_indiv)))
            elseif g == 2
                global timings[index] = @elapsed reg(base, @formula(y ~ x1 + x2 + fe(id_indiv) + fe(id_firm)))
            else
                global timings[index] = @elapsed reg(base, @formula(y ~ x1 + x2 + fe(id_indiv) + fe(id_firm) + fe(id_year)))
            end
            global index = index + 1
        end
    end
end

# saving
open("DATA/julia_bench_diff.txt", "w") do io
           writedlm(io, timings)
       end
