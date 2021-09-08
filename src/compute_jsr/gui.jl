Cursor = matplotlib.widgets.Cursor
Button = matplotlib.widgets.Button
CheckButtons = matplotlib.widgets.CheckButtons
Slider = matplotlib.widgets.Slider

mutable struct the_gui_compute_jsr{N}
    A_list::NTuple{N,SMatrix{2,2,Float64}}
    γ::Float64
    trace::Float64
    ηm::Float64
    ηa::Float64
    P_opt::SMatrix{2,2,Float64}
    is_P_opt_set::Bool
    is_singular::Bool
    index::Int
    fig
    ax_points
    ax_ellipse
    ax_constr
    plt_points_x
    plt_points_ηm
    plt_points_ηa
    plt_points_y
    plt_points_indexes
    plt_ellipse
    plt_ellipse_γ
    plt_ellipse_x
    plt_ellipse_y
    plt_constr_circle
    plt_constr_ellipse
    plt_constr_lines
    x_list::Vector{SVector{2,Float64}}
    y_list::Vector{SVector{2,Float64}}
    ηm_list::Vector{Float64}
    ηa_list::Vector{Float64}
    cursor_points
    button_compute
    button_clear
    button_infty
    slider_γ
    slider_tr
    slider_ηm
    slider_ηa
    slider_Ai
    circle::Tuple{Vector{Float64},Vector{Float64}}
end

function new_gui_compute_jsr(A_list::NTuple{N,SMatrix{2,2,Float64}};
        γ_max=5.0, trace_max=50.0, γ=1.0, trace=20.0, np=200) where N
    @assert N > 0
    fig = PyPlot.figure(figsize = (10.4, 9.8))
    ax_points = fig.add_axes((0.05, 0.57, 0.4, 0.4), aspect = "equal")
    ax_ellipse = fig.add_axes((0.55, 0.57, 0.4, 0.4), aspect = "equal")
    ax_constr = fig.add_axes((0.05, 0.05, 0.4, 0.4), aspect = "equal")

    # Set useblit=True on most backends for enhanced performance.
    cursor_points = Cursor(ax_points, useblit = true, color = "red", linewidth = 1.5)
    plt_points_x, = ax_points.plot((), ())
    plt_points_ηm = ax_points.scatter((), (), (), marker = "x")
    plt_points_ηa = ax_points.scatter((), (), (), marker = "o")
    plt_points_y, = ax_points.plot((), ())
    plt_points_indexes = []

    plt_ellipse, = ax_ellipse.plot((), ())
    plt_ellipse_γ, = ax_ellipse.plot((), ())
    plt_ellipse_x, = ax_ellipse.plot((), ())
    plt_ellipse_y, = ax_ellipse.plot((), ())

    plt_constr_circle, = ax_constr.plot((), ())
    plt_constr_ellipse, = ax_constr.plot((), ())
    plt_constr_lines = []

    ax_button_compute = fig.add_axes((0.76, 0.33, 0.13, 0.075))
    button_compute = Button(ax_button_compute, "Compute cone")
    ax_button_clear = fig.add_axes((0.63, 0.33, 0.1, 0.075))
    button_clear = Button(ax_button_clear, "Clear")
    ax_button_infty = fig.add_axes((0.51, 0.3225, 0.09, 0.09))
    button_infty = CheckButtons(ax_button_infty, ("trace inf",), (false,))

    ax_slider_γ = fig.add_axes((0.06, 0.48, 0.88, 0.03))
    slider_γ = Slider(ax_slider_γ, "γ", valmin=0.0, valmax=γ_max, valinit=γ)
    ax_slider_tr = fig.add_axes((0.5, 0.23, 0.4, 0.03))
    slider_tr = Slider(ax_slider_tr, "tr", valmin=0.0, valmax=trace_max, valinit=trace)
    ax_slider_ηm = fig.add_axes((0.5, 0.19, 0.4, 0.03))
    slider_ηm = Slider(ax_slider_ηm, "ηm", valmin=0.0, valmax=1.0, valinit=0.0)
    ax_slider_ηa = fig.add_axes((0.5, 0.15, 0.4, 0.03))
    slider_ηa = Slider(ax_slider_ηa, "ηa", valmin=0.0, valmax=1.0, valinit=0.0)
    ax_slider_Ai = fig.add_axes((0.5, 0.11, 0.4, 0.03))
    slider_Ai = Slider(ax_slider_Ai, L"A_i", 1, N, valstep=1, valinit=1)

    the_t_circle = range(0.0, 2.0*π, length = np)

    my_gui = the_gui_compute_jsr(
        A_list, γ, trace, 0.0, 0.0, zero(SMatrix{2,2}), false, false, 1,
        fig, ax_points, ax_ellipse, ax_constr,
        plt_points_x, plt_points_ηm, plt_points_ηa, plt_points_y, plt_points_indexes,
        plt_ellipse, plt_ellipse_γ, plt_ellipse_x, plt_ellipse_y,
        plt_constr_circle, plt_constr_ellipse, plt_constr_lines,
        SVector{2,Float64}[], SVector{2,Float64}[], Float64[], Float64[],
        cursor_points, button_compute, button_clear, button_infty,
        slider_γ, slider_tr, slider_ηm, slider_ηa, slider_Ai,
        (cos.(the_t_circle), sin.(the_t_circle)))

    initialize(my_gui)

    return my_gui
end

function initialize(MG::the_gui_compute_jsr)
    MG.ax_points.plot(MG.circle[1], MG.circle[2], "tab:blue")
    MG.ax_points.set_xlim(-1.5, 1.5)
    MG.ax_points.set_ylim(-1.5, 1.5)
    MG.ax_ellipse.set_xlim(-1.5, 1.5)
    MG.ax_ellipse.set_ylim(-1.5, 1.5)

    MG.plt_points_x.set_ls("none")
    MG.plt_points_x.set_marker("o")
    MG.plt_points_x.set_c("blue")
    MG.plt_points_ηm.set_ec("orange")
    MG.plt_points_ηm.set_fc("none")
    MG.plt_points_ηm.set_lw(1.5)
    MG.plt_points_ηm.set_zorder(100)
    MG.plt_points_ηa.set_ec("green")
    MG.plt_points_ηa.set_fc("none")
    MG.plt_points_ηa.set_lw(1.5)
    MG.plt_points_ηa.set_zorder(100)
    MG.plt_points_y.set_ls("none")
    MG.plt_points_y.set_marker("o")
    MG.plt_points_y.set_c("red")

    MG.plt_ellipse.set_lw("2.0")
    MG.plt_ellipse.set_c("black")
    MG.plt_ellipse_γ.set_ls("--")
    MG.plt_ellipse_γ.set_c("black")
    MG.plt_ellipse_x.set_ls("none")
    MG.plt_ellipse_x.set_marker("o")
    MG.plt_ellipse_x.set_c("blue")
    MG.plt_ellipse_y.set_ls("none")
    MG.plt_ellipse_y.set_marker("o")
    MG.plt_ellipse_y.set_c("red")

    MG.plt_constr_ellipse.set_ls("none")
    MG.plt_constr_ellipse.set_marker("o")
    MG.plt_constr_ellipse.set_c("red")
    draw_constr_circle(MG)
    set_constr_lims(MG)

    MG.button_compute.on_clicked() do event
        on_press_button_compute(MG)
    end

    MG.button_clear.on_clicked() do event
        on_press_button_clear(MG)
    end

    MG.button_infty.on_clicked() do event
        MG.is_singular = MG.button_infty.get_status()[1]
        draw_constr_circle(MG)
        if MG.is_P_opt_set
            draw_constr_ellipse(MG)
        end
        draw_constr_lines(MG)
        set_constr_lims(MG)
    end

    MG.slider_γ.on_changed() do val
        MG.γ = val
        draw_constr_lines(MG)
    end

    MG.slider_tr.on_changed() do val
        MG.trace = val
        if !MG.is_singular
            draw_constr_circle(MG)
            if MG.is_P_opt_set
                draw_constr_ellipse(MG)
            end
            draw_constr_lines(MG)
            set_constr_lims(MG)
        end
    end

    MG.slider_ηm.on_changed() do val
        MG.ηm = val
    end

    MG.slider_ηa.on_changed() do val
        MG.ηa = val
    end

    MG.slider_Ai.on_changed() do val
        MG.index = val
    end

    MG.fig.canvas.mpl_connect("button_press_event", event -> begin
        @printf("%s click: button=%d, x=%d, y=%d\n",
            event.dblclick ? "double" : "single", event.button, event.x, event.y)

        if event.inaxes == MG.ax_points
            xydata = (event.xdata, event.ydata)
            on_press_points(xydata, MG)
            return
        end
        if event.inaxes == MG.button_compute.ax ||
                event.inaxes == MG.button_clear.ax ||
                event.inaxes == MG.slider_γ.ax ||
                event.inaxes == MG.slider_tr.ax ||
                event.inaxes == MG.slider_ηm.ax ||
                event.inaxes == MG.slider_ηa.ax
            @printf("pressed button/slider\n")
            return
        end

        @printf("void place: %s\n", event.inaxes)
    end)
end

function on_press_points(xydata, MG)
    x = SVector(xydata)
    y = MG.A_list[MG.index]*x
    push!(MG.x_list, x)
    push!(MG.ηm_list, MG.ηm)
    push!(MG.ηa_list, MG.ηa)
    push!(MG.y_list, y)
    MG.plt_points_x.set_xdata(map(p -> p[1], MG.x_list))
    MG.plt_points_x.set_ydata(map(p -> p[2], MG.x_list))
    MG.plt_points_ηm.set_offsets(MG.x_list)
    MG.plt_points_ηm.set_sizes(MG.ηm_list.*1000)
    MG.plt_points_ηa.set_offsets(MG.x_list)
    MG.plt_points_ηa.set_sizes(MG.ηa_list.*1000)
    MG.plt_points_y.set_xdata(map(p -> p[1], MG.y_list))
    MG.plt_points_y.set_ydata(map(p -> p[2], MG.y_list))
    push!(MG.plt_points_indexes, MG.ax_points.text(x[1], x[2], MG.index))
    the_trace = MG.is_singular ? 1.0 : MG.trace
    p, q = points_2_plane(x, y, MG.γ, the_trace, MG.ηm, MG.ηa)
    PLT, = MG.ax_constr.plot((p[1], q[1]), (p[2], q[2]))
    PLT.set_c("black")
    push!(MG.plt_constr_lines, PLT)
end

function on_press_button_compute(MG)
    the_max_trace = MG.is_singular ? nothing : 1e5
    γ_opt, P_opt = jsr_quadratic(
        MG.x_list, MG.y_list, MG.ηm_list, MG.ηa_list, 0.0, 1e2;
        max_trace=the_max_trace)
    MG.P_opt = P_opt
    MG.is_P_opt_set = true
    @printf("Trace: %f\n", tr(P_opt))
    n_Q, x1_ellipse_list, x2_ellipse_list = matrix_2_ellipse(P_opt, 500)
    MG.plt_ellipse.set_xdata(x1_ellipse_list)
    MG.plt_ellipse.set_ydata(x2_ellipse_list)
    MG.plt_ellipse_γ.set_xdata(γ_opt.*x1_ellipse_list)
    MG.plt_ellipse_γ.set_ydata(γ_opt.*x2_ellipse_list)

    xplist = SVector{2,Float64}[]
    yplist = SVector{2,Float64}[]

    for (x, y) in zip(MG.x_list, MG.y_list)
        n_x = sqrt(x'*P_opt*x)*n_Q
        push!(xplist, x/n_x)
        push!(yplist, y/n_x)
    end

    MG.plt_ellipse_x.set_xdata(map(p -> p[1], xplist))
    MG.plt_ellipse_x.set_ydata(map(p -> p[2], xplist))
    MG.plt_ellipse_y.set_xdata(map(p -> p[1], yplist))
    MG.plt_ellipse_y.set_ydata(map(p -> p[2], yplist))
    draw_constr_ellipse(MG)
end

function on_press_button_clear(MG)
    empty!(MG.x_list)
    empty!(MG.y_list)
    empty!(MG.ηm_list)
    empty!(MG.ηa_list)
    MG.is_P_opt_set = false
    MG.plt_points_x.set_xdata(())
    MG.plt_points_x.set_ydata(())
    MG.plt_points_ηm.set_offsets(())
    MG.plt_points_ηm.set_sizes(())
    MG.plt_points_ηa.set_offsets(())
    MG.plt_points_ηa.set_sizes(())
    MG.plt_points_y.set_xdata(())
    MG.plt_points_y.set_ydata(())
    MG.ax_points.set_xlim(-1.5, 1.5)
    MG.ax_points.set_ylim(-1.5, 1.5)
    MG.ax_ellipse.set_xlim(-1.5, 1.5)
    MG.ax_ellipse.set_ylim(-1.5, 1.5)
    set_constr_lims(MG)
    for slider in (MG.slider_γ, MG.slider_tr, MG.slider_ηm, MG.slider_ηa, MG.slider_Ai)
        slider.reset()
    end
    for PLT in MG.plt_points_indexes
        PLT.remove()
    end
    empty!(MG.plt_points_indexes)
    for PLT in MG.plt_constr_lines
        PLT.remove()
    end
    empty!(MG.plt_constr_lines)
end

function draw_constr_lines(MG)
    the_trace = MG.is_singular ? 1.0 : MG.trace
    for (x, y, ηm, ηa, PLT) in zip(MG.x_list, MG.y_list,
            MG.ηm_list, MG.ηa_list, MG.plt_constr_lines)
        p, q = points_2_plane(x, y, MG.γ, the_trace, ηm, ηa)
        PLT.set_xdata((p[1], q[1]))
        PLT.set_ydata((p[2], q[2]))
    end
end

function draw_constr_circle(MG)
    rad = MG.is_singular ? 0.5 : max(0.0, MG.trace/2.0 - 1.0)
    MG.plt_constr_circle.set_xdata(rad.*MG.circle[1])
    MG.plt_constr_circle.set_ydata(rad.*MG.circle[2])
end

function draw_constr_ellipse(MG)
    the_trace = MG.is_singular ? 1.0 : MG.trace
    P_opt_n = the_trace*MG.P_opt./tr(MG.P_opt)
    t = (P_opt_n[1, 1] - P_opt_n[2, 2])/2
    u = (P_opt_n[1, 2] + P_opt_n[2, 1])/2
    MG.plt_constr_ellipse.set_xdata((t,))
    MG.plt_constr_ellipse.set_ydata((u,))
end

function set_constr_lims(MG)
    clims = MG.is_singular ? 0.5 : max(1.0, MG.trace/2.0 - 1.0)
    MG.ax_constr.set_xlim(-1.05*clims, 1.05*clims)
    MG.ax_constr.set_ylim(-1.05*clims, 1.05*clims)
end

function matrix_2_ellipse(P, np)
    EV = eigen(Symmetric(Matrix(P)))
    EVq = Eigen(sqrt.(EV.values), EV.vectors)
    Q = SMatrix{2,2}(inv(EVq))
    n_Q = opnorm(Q)
    Q = Q/n_Q
    t_list = range(0.0, 2.0*π, length = np)
    x_list = map(t -> Q*SVector(cos(t), sin(t)), t_list)
    return n_Q, map(p -> p[1], x_list), map(p -> p[2], x_list)
end

function points_2_plane(x, y, γ, C, ηm, ηa)
    δ = γ^2*(1.0 + ηm)
    v = SVector(y[1]^2 - y[2]^2 - δ*(x[1]^2 - x[2]^2), 2*y[1]*y[2] - 2*δ*x[1]*x[2])
    nv = norm(v)
    b = -0.5*C*(y[1]^2 + y[2]^2 - δ*(x[1]^2 + x[2]^2)) + δ*ηa
    z = v*(b/nv^2)
    w = SVector(v[2], -v[1])/nv
    rad = max(0.0, C/2.0 - 1.0)
    length = max(1.0, rad)
    return z + length*w, z - length*w
end