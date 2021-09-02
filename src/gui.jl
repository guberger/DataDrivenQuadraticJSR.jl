Cursor = matplotlib.widgets.Cursor
Button = matplotlib.widgets.Button
Slider = matplotlib.widgets.Slider

function new_gui(A)
    fig = PyPlot.figure(figsize = (10.4, 9.8))
    ax_points = fig.add_axes((0.05, 0.57, 0.4, 0.4), aspect = "equal")
    ax_ellipse = fig.add_axes((0.55, 0.57, 0.4, 0.4), aspect = "equal")
    ax_constr = fig.add_axes((0.05, 0.05, 0.4, 0.4), aspect = "equal")

    the_γ = [1.0]
    the_tr = [20.0]
    the_ηm = [0.0]
    the_ηa = [0.0]
    the_P_opt = [SMatrix{2,2}(NaN, NaN, NaN, NaN)]

    np = 200
    the_t_circle = range(0.0, 2.0*π, length = np)
    ax_points.plot(cos.(the_t_circle), sin.(the_t_circle))
    rad = max(0.0, the_tr[1]/2.0 - 1.0)
    plt_constr_circle, = ax_constr.plot(
        rad.*cos.(the_t_circle), rad.*sin.(the_t_circle))
    ax_points.set_xlim(-1.5, 1.5)
    ax_points.set_ylim(-1.5, 1.5)
    ax_ellipse.set_xlim(-1.5, 1.5)
    ax_ellipse.set_ylim(-1.5, 1.5)
    clims = max(1.0, rad)
    ax_constr.set_xlim(-1.05*clims, 1.05*clims)
    ax_constr.set_ylim(-1.05*clims, 1.05*clims)
    # Set useblit=True on most backends for enhanced performance.
    cursor_points = Cursor(ax_points, useblit = true, color = "red",
        linewidth = 1.5)
    plt_points_x, = ax_points.plot((), ())
    plt_points_x.set_ls("none")
    plt_points_x.set_marker("o")
    plt_points_x.set_c("blue")
    plt_points_ηm = ax_points.scatter((), (), (), marker = "x")
    plt_points_ηm.set_ec("orange")
    plt_points_ηm.set_fc("none")
    plt_points_ηm.set_lw(1.5)
    plt_points_ηm.set_zorder(100)
    plt_points_ηa = ax_points.scatter((), (), (), marker = "o")
    plt_points_ηa.set_ec("green")
    plt_points_ηa.set_fc("none")
    plt_points_ηa.set_lw(1.5)
    plt_points_ηa.set_zorder(100)
    plt_points_y, = ax_points.plot((), ())
    plt_points_y.set_ls("none")
    plt_points_y.set_marker("o")
    plt_points_y.set_c("red")

    plt_ellipse, = ax_ellipse.plot((), ())
    plt_ellipse.set_lw("2.0")
    plt_ellipse.set_c("black")
    plt_ellipse_γ, = ax_ellipse.plot((), ())
    plt_ellipse_γ.set_ls("--")
    plt_ellipse_γ.set_c("black")
    plt_ellipse_x, = ax_ellipse.plot((), ())
    plt_ellipse_x.set_ls("none")
    plt_ellipse_x.set_marker("o")
    plt_ellipse_x.set_c("blue")
    plt_ellipse_y, = ax_ellipse.plot((), ())
    plt_ellipse_y.set_ls("none")
    plt_ellipse_y.set_marker("o")
    plt_ellipse_y.set_c("red")

    plt_constr_ellipse, = ax_constr.plot((), ())
    plt_constr_ellipse.set_ls("none")
    plt_constr_ellipse.set_marker("o")
    plt_constr_ellipse.set_c("red")
    plt_constr_list = []

    the_x_list = SVector{2,Float64}[]
    the_y_list = SVector{2,Float64}[]
    the_ηm_list = Float64[]
    the_ηa_list = Float64[]

    ax_button_compute = fig.add_axes((0.76, 0.33, 0.13, 0.075))
    button_compute = Button(ax_button_compute, "Compute cone")
    ax_button_clear = fig.add_axes((0.63, 0.33, 0.1, 0.075))
    button_clear = Button(ax_button_clear, "Clear")

    button_compute.on_clicked() do event
        on_press_button_compute(the_x_list, the_y_list, the_ηm_list, the_ηa_list,
            plt_ellipse, plt_ellipse_γ, plt_ellipse_x, plt_ellipse_y,
            plt_constr_ellipse, the_tr, the_P_opt)
    end

    button_clear.on_clicked() do event
        on_press_button_clear(the_x_list, the_y_list, the_ηm_list, the_ηa_list,
            ax_points, ax_ellipse, ax_constr,
            plt_points_x, plt_points_ηm, plt_points_ηa, plt_points_y,
            plt_constr_list,
            (slider_γ, slider_tr, slider_ηm, slider_ηa), the_tr, the_P_opt)
    end

    ax_slider_γ = fig.add_axes((0.06, 0.48, 0.88, 0.03))
    slider_γ = Slider(ax_slider_γ, "γ", valmin = 0.0, valmax = 5.0,
        valinit = the_γ[1])
    ax_slider_tr = fig.add_axes((0.5, 0.23, 0.4, 0.03))
    slider_tr = Slider(ax_slider_tr, "tr", valmin = 0.0, valmax = 50.0,
        valinit = the_tr[1])
    ax_slider_ηm = fig.add_axes((0.5, 0.19, 0.4, 0.03))
    slider_ηm = Slider(ax_slider_ηm, "ηm", valmin = 0.0, valmax = 1.0,
        valinit = the_ηm[1])
    ax_slider_ηa = fig.add_axes((0.5, 0.15, 0.4, 0.03))
    slider_ηa = Slider(ax_slider_ηa, "ηa", valmin = 0.0, valmax = 1.0,
        valinit = the_ηa[1])

    slider_γ.on_changed() do val
        the_γ[1] = val
        draw_constraints(the_x_list, the_y_list, the_ηm_list, the_ηa_list,
            plt_constr_list, the_γ[1], the_tr[1])
    end

    slider_tr.on_changed() do val
        the_tr[1] = val
        rad = max(0.0, the_tr[1]/2.0 - 1.0)
        plt_constr_circle.set_xdata(rad.*cos.(the_t_circle))
        plt_constr_circle.set_ydata(rad.*sin.(the_t_circle))
        draw_constraints(the_x_list, the_y_list, the_ηm_list, the_ηa_list,
            plt_constr_list, the_γ[1], the_tr[1])
        clims = max(1.0, rad)
        ax_constr.set_xlim(-1.05*clims, 1.05*clims)
        ax_constr.set_ylim(-1.05*clims, 1.05*clims)
        if the_P_opt[1] == SMatrix{2,2}(NaN, NaN, NaN, NaN)
            return
        end
        P_opt = the_P_opt[1]
        P_opt_n = the_tr[1]*P_opt./tr(P_opt)
        t = (P_opt_n[1, 1] - P_opt_n[2, 2])/2
        u = (P_opt_n[1, 2] + P_opt_n[2, 1])/2
        plt_constr_ellipse.set_xdata((t,))
        plt_constr_ellipse.set_ydata((u,))
    end

    slider_ηm.on_changed() do val
        the_ηm[1] = val
    end

    slider_ηa.on_changed() do val
        the_ηa[1] = val
    end

    fig.canvas.mpl_connect("button_press_event", event -> begin
        @printf("%s click: button=%d, x=%d, y=%d\n",
            event.dblclick ? "double" : "single", event.button, event.x, event.y)

        if event.inaxes == ax_points
            xydata = (event.xdata, event.ydata)
            on_press_points(xydata, A,
                the_x_list, the_y_list, the_ηm_list, the_ηa_list,
                ax_constr,
                plt_points_x, plt_points_ηm, plt_points_ηa, plt_points_y,
                plt_constr_list,
                the_γ, the_tr, the_ηm, the_ηa)
            return
        end
        if event.inaxes == ax_button_compute || event.inaxes == ax_button_clear
            @printf("pressed button\n")
            return
        end

        @printf("void place: %s\n", event.inaxes)
    end)

    return cursor_points, button_compute, button_clear,
        slider_γ, slider_tr, slider_ηm, slider_ηa,
        the_P_opt
end

function on_press_points(xydata, A, xlist, ylist, ηmlist, ηalist,
    axc, pltx, pltηm, pltηa, plty, pltclist, the_γ, the_tr, the_ηm, the_ηa)
    x = SVector(xydata)
    y = A*x
    push!(xlist, x)
    push!(ηmlist, the_ηm[1])
    push!(ηalist, the_ηa[1])
    push!(ylist, y)
    pltx.set_xdata(map(p -> p[1], xlist))
    pltx.set_ydata(map(p -> p[2], xlist))
    pltηm.set_offsets(xlist)
    pltηm.set_sizes(ηmlist.*1000)
    pltηa.set_offsets(xlist)
    pltηa.set_sizes(ηalist.*1000)
    plty.set_xdata(map(p -> p[1], ylist))
    plty.set_ydata(map(p -> p[2], ylist))
    p, q = points_2_plane(x, y, the_γ[1], the_tr[1], the_ηm[1], the_ηa[1])
    PLT, = axc.plot((p[1], q[1]), (p[2], q[2]))
    PLT.set_c("black")
    push!(pltclist, PLT)
end

function on_press_button_compute(xlist, ylist, ηmlist, ηalist,
        plte, plteγ, pltex, pltey, pltce, the_tr, the_P_opt)
    γ_opt, P_opt = jsr_quadratic(xlist, ylist, ηmlist, ηalist, 0.0, 1e2)
    the_P_opt[1] = P_opt
    trace = tr(P_opt)
    n_Q, x1_ellipse_list, x2_ellipse_list = matrix_2_ellipse(P_opt, 500)
    plte.set_xdata(x1_ellipse_list)
    plte.set_ydata(x2_ellipse_list)
    plteγ.set_xdata(γ_opt.*x1_ellipse_list)
    plteγ.set_ydata(γ_opt.*x2_ellipse_list)

    xplist = SVector{2,Float64}[]
    yplist = SVector{2,Float64}[]

    for (x, y) in zip(xlist, ylist)
        n_x = sqrt(x'*P_opt*x)*n_Q
        push!(xplist, x/n_x)
        push!(yplist, y/n_x)
    end

    pltex.set_xdata(map(p -> p[1], xplist))
    pltex.set_ydata(map(p -> p[2], xplist))
    pltey.set_xdata(map(p -> p[1], yplist))
    pltey.set_ydata(map(p -> p[2], yplist))

    P_opt_n = the_tr[1]*P_opt./tr(P_opt)
    t = (P_opt_n[1, 1] - P_opt_n[2, 2])/2
    u = (P_opt_n[1, 2] + P_opt_n[2, 1])/2
    pltce.set_xdata((t,))
    pltce.set_ydata((u,))
    @printf("Trace: %f\n", trace)
end

function on_press_button_clear(xlist, ylist, ηmlist, ηalist,
    axp, axe, axc, pltx, pltηm, pltηa, plty, pltclist, sliders, the_tr, the_P_opt)
    empty!(xlist)
    empty!(ylist)
    empty!(ηmlist)
    empty!(ηalist)
    the_P_opt[1] = SMatrix{2,2}(NaN, NaN, NaN, NaN)
    pltx.set_xdata(())
    pltx.set_ydata(())
    pltηm.set_offsets(())
    pltηm.set_sizes(())
    pltηa.set_offsets(())
    pltηa.set_sizes(())
    plty.set_xdata(())
    plty.set_ydata(())
    axp.set_xlim(-1.5, 1.5)
    axp.set_ylim(-1.5, 1.5)
    axe.set_xlim(-1.5, 1.5)
    axe.set_ylim(-1.5, 1.5)
    rad = max(0.0, the_tr[1]/2.0 - 1.0)
    clims = max(1.0, rad)
    axc.set_xlim(-1.05*clims, 1.05*clims)
    axc.set_ylim(-1.05*clims, 1.05*clims)
    for slider in sliders
        slider.reset()
    end
    for PLT in pltclist
        PLT.remove()
    end
    empty!(pltclist)
end

function draw_constraints(xlist, ylist, ηmlist, ηalist, pltclist, γ, C)
    for (x, y, ηm, ηa, PLT) in zip(xlist, ylist, ηmlist, ηalist, pltclist)
        p, q = points_2_plane(x, y, γ, C, ηm, ηa)
        PLT.set_xdata((p[1], q[1]))
        PLT.set_ydata((p[2], q[2]))
    end
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
    clims = max(1.0, rad)
    return z + clims*w, z - clims*w
end