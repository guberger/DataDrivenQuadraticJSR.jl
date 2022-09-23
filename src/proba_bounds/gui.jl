Slider = matplotlib.widgets.Slider

mutable struct the_gui_proba_bounds
    ax_plots
    plot1
    plot2
    plot3
    slider_β
    slider_N
    slider_ζ
    logβ::Float64
    N::Int
    ζ::Int
    tol::Float64
end

function new_gui_proba_bounds(; logβmin::Float64=-10.0, Nmax::Int=500, ζ = 10, tol = 1e-6)
    @assert Nmax > 1 && logβmin < 0
    logβ = logβmin/2
    N = Nmax
    ζ = max(0, min(ζ, Nmax))

    fig = PyPlot.figure(figsize = (14.4, 7.0))
    ax_plots = fig.add_axes((0.1, 0.07, 0.6, 0.5))
    plot1, = ax_plots.plot((), ())
    plot2, = ax_plots.plot((), ())
    plot3, = ax_plots.plot((), ())

    ax_slider_β = fig.add_axes((0.08, 0.88, 0.84, 0.03))
    slider_β = Slider(ax_slider_β, L"log_{10}(\beta)", logβmin, 0.0, valinit=logβ)
    ax_slider_β = fig.add_axes((0.08, 0.78, 0.84, 0.03))
    slider_N = Slider(ax_slider_β, L"N", 1, Nmax, valstep=1, valinit=N)
    ax_slider_ζ = fig.add_axes((0.08, 0.68, 0.84, 0.03))
    slider_ζ = Slider(ax_slider_ζ, L"\zeta", 0, Nmax, valstep=1, valinit=ζ)

    my_gui = the_gui_proba_bounds(ax_plots, plot1, plot2, plot3,
        slider_β, slider_N, slider_ζ,
        logβ, N, ζ, tol)

    initialize(my_gui)
    
    return my_gui
end

function initialize(MG::the_gui_proba_bounds)
    MG.plot1.set_marker("o")
    MG.plot1.set_c("tab:blue")
    MG.plot2.set_marker("x")
    MG.plot2.set_c("tab:green")
    MG.plot3.set_marker("d")
    MG.plot3.set_c("tab:orange")

    LH = (
        matplotlib.lines.Line2D([0], [0], marker = "o", c = "tab:blue",
            label = "Uniform"),
        matplotlib.lines.Line2D([0], [0], marker = "x", c = "tab:green",
            label = "Varying 2018"),
        matplotlib.lines.Line2D([0], [0], marker = "d", c = "tab:orange",
            label = "Varying 2021"))
    MG.ax_plots.legend(handles = LH, loc = "upper left", bbox_to_anchor = (1.01, 1.0))

    draw_plots(MG)

    MG.slider_β.on_changed() do val
        MG.logβ = val
        draw_plots(MG)
    end
    MG.slider_N.on_changed() do val
        MG.N = val
        MG.slider_ζ.valmax = val
        MG.slider_ζ.ax.set_xlim((MG.slider_ζ.valmin, val))
        ζ = min(MG.ζ, val)
        MG.slider_ζ.set_val(ζ)
        if ζ == MG.ζ
            draw_plots(MG)
        end
    end
    MG.slider_ζ.on_changed() do val
        MG.ζ = val
        draw_plots(MG)
    end
end

function draw_plots(MG)
    k_list = 0:MG.ζ
    ~, ϵ1 = uniform_bound(MG.N, MG.ζ, exp10(MG.logβ), tol=MG.tol)
    MG.plot1.set_xdata(k_list)
    MG.plot1.set_ydata(fill(ϵ1, length(k_list)))
    ϵ2_list = map(k -> varying_bound_1(MG.N, k, exp10(MG.logβ), tol=MG.tol)[2], k_list)
    MG.plot2.set_xdata(k_list)
    MG.plot2.set_ydata(ϵ2_list)
    ϵ3_list = map(k -> varying_bound_2(MG.N, k, exp10(MG.logβ), tol=MG.tol)[2], k_list)
    MG.plot3.set_xdata(k_list)
    MG.plot3.set_ydata(ϵ3_list)
    xmax = max(MG.ζ, 1)
    MG.ax_plots.set_xlim((0, xmax))
    ymax = max(ϵ1, maximum(ϵ2_list), maximum(ϵ3_list))*1.2
    MG.ax_plots.set_ylim((0, ymax))
end