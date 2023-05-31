
using Pkg; 
Pkg.add("GaussianRandomFields"); Pkg.add("Plots"); Pkg.add("MAT");
using GaussianRandomFields;
using Plots;
using MAT 

# 3D scenario
n_dimension = 3 
nx = 100 
ny = 100
nz = 11
nsamples = 101

# cov = CovarianceFunction(n_dimension, Exponential(0.1, σ=0.5)) # σ=0.6
cov = CovarianceFunction(n_dimension, Exponential(25, σ=3.0, p=4.0)) # Exponential(correlation length, std, p-norm)
xpts = range(0, stop=1, length=nx)
ypts = range(0, stop=1, length=ny)
zpts = range(0, stop=1, length=nz)
# using a Karhunen-Loève (KL) expansion with n terms, n=20 in Dongxiao Papper
grf_perm = GaussianRandomField(3, cov, KarhunenLoeve(20), xpts, ypts, zpts) # μ=3.0, 
heatmap(exp(1) .^(sample(grf_perm)))

perm = zeros(Float64, nx, ny, nz, nsamples)
for i=1:nsamples
    for j=1:10000
        a = 5 .^(sample(grf_perm))
        if 0 < minimum(a) < 20 && 50 < maximum(a) < 150
            break
        end
    end
    perm[:, :, :, i] = a
end

matwrite(
    "PERMALL_KL_n_100_100_11_101.mat",
    Dict("PERMALL" => perm)
)

