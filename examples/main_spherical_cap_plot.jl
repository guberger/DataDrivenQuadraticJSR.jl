using PyPlot
using Distributions
_cols = repeat(matplotlib.rcParams["axes.prop_cycle"].by_key()["color"], 10, 1)
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

fig = figure(figsize = (18.0, 8.5))
ax = fig.add_subplot()

η_list = range(0, 1/2, length=300)
dim_list = 2:8

for (i, n) in enumerate(dim_list)
    s_list = map(η -> DataDrivenQuadraticJSR.area_2sided_cap_inv(2*η, n), η_list)
    ax.plot(η_list, s_list, c=_cols[i], ls="-", lw=2.5)
    # Ap_list = map(s -> DataDrivenQuadraticJSR.area_2sided_cap(s, n), s_list)
    # ax.plot(Ap_list, s_list, c="k", ls="--", lw=2.5)
end

ax.set_xlabel(L"$\eta$")
ax.set_ylabel(L"$\bar{s}$")

LH = [matplotlib.lines.Line2D([0], [0], c=_cols[i], ls="-", lw=4,
        label = "\$n=$(dim_list[i])\$") for i in eachindex(dim_list)]
ax.legend(handles = LH)

fig.savefig(string("./figures/fig_spherical_cap.png"), dpi=200,
    transparent = false, bbox_inches = "tight")

print("")
sleep(0.1)