using Observables
using GLMakie

include("testpolymeranalysis.jl")


function replot(raw_data_file)

    x, y, v_x, v_y, p_x, p_y, pol_id, numb_frames, numb_particles, dt, Tsave = get_info(raw_data_file)

    x_o = Observable(x[1])
    y_o = Observable(y[1])

    p_x_o = Observable(p_x[1])
    p_y_o = Observable(p_y[1])

    t0 = Observable(0.0)
    
    fig = Figure([1,1])
    ax = Axis(f[1, 1], xlabel = "x", ylabel="y", title = @lift("t = $($tO)"))


    for i in 2:numb_frames
        x_o[] = x[i]
        y_o[] = y[i]
        p_x_o[] = p_x[i]
        p_y_o[] = p_y[i]
        t0[] = i*dt*Tsave

        x = x_o
        y = y_o
        p = [p_x_o, p_y_o]
        pmin = [-p_x_o, -p_y_o]

        scatter!(ax,x,y, color=pol_id, markersize =2,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1, colormap=:rainbow1)
        arrows!(ax, x, p , color=c,  colormap=:hsv,colorrange=(-pi,pi))
        arrows!(ax, x, pmin , color=c,  colormap=:hsv,colorrange=(-pi,pi))
        display(fig)
    end


end

file = load_file("C:\\Users\\gabri\\Documents\\Travail-Etude\\Master's Theoretical Physics Leiden\\Research Project\\Data\\sim_data\\p_0.13\\raw_data.h5")

replot(file)