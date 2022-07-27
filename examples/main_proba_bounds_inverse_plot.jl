using PyPlot
_cols = repeat(matplotlib.rcParams["axes.prop_cycle"].by_key()["color"], 10, 1)
_styles = repeat(["-", "--", ".", "-."], 10, 1)
matplotlib.rc("legend", fontsize = 30)
matplotlib.rc("axes", labelsize = 30)
matplotlib.rc("xtick", labelsize = 30)
matplotlib.rc("ytick", labelsize = 30)
matplotlib.rc("text", usetex = true)
matplotlib.rc("text.latex", preamble = "\\usepackage{amsmath,amssymb}")

include(joinpath(@__DIR__, "../src/DataDrivenQuadraticJSR.jl"))
# using Main.DataDrivenQuadraticJSR

println("\nNew test")
sleep(0.1)

β_list = (0.01, 0.001)
ϵ_list = 0.5 .^ (1:10)
n_list = 1:4

fig = figure(figsize = (18.0, 8.5))
ax = fig.add_subplot()
ax.set_xscale("log")
ax.set_yscale("log")

for (i, β) in enumerate(β_list)
    for (j, n) in enumerate(n_list)
        # local k = n^2 + n - 1
        local k = Int(n*(n + 1)/2)
        local N_list = zeros(Int, length(ϵ_list))
        for (r, ϵ) in enumerate(ϵ_list)
            N_list[r] = DataDrivenQuadraticJSR.varying_bound_2_inv(β, ϵ, k)
        end
        ax.plot(N_list, ϵ_list, c=_cols[j], ls=_styles[i], lw=2.5)
    end
    # local N_list = zeros(Int, length(ϵ_list))
    # for (r, ϵ) in enumerate(ϵ_list)
    #     N_list[r] = DataDrivenQuadraticJSR.single_monotone_bound_inv(β, ϵ)
    # end
    # ax.plot(N_list, ϵ_list, c=_cols[length(n_list) + 1], ls=_styles[i], lw=2.5)
end

ax.set_xlabel(L"$N$")
# ax.set_ylabel(L"$\bar\epsilon(n^2+n-1)$")
ax.set_ylabel(L"$\bar\epsilon(\frac{n(n+1)}2)$")

LH = vcat(
    [matplotlib.lines.Line2D([0], [0], c="k", ls=_styles[i], lw=4,
        label = "\$\\beta=$(β_list[i])\$") for i in eachindex(β_list)],
    [matplotlib.patches.Patch(fc=_cols[j], ec="none",
        label = "\$n=$(n_list[j])\$") for j in eachindex(n_list)])
ax.legend(handles = LH)

fig.savefig(string("./figures/fig_proba_bounds.png"), dpi=200,
    transparent = false, bbox_inches = "tight")

print("")
sleep(0.1)