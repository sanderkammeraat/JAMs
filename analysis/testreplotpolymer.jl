using Observables
using GLMakie

include("testpolymeranalysis.jl")

GLMakie.activate()

function replot(raw_data_file)

    x, y, v_x, v_y, numb_frames, numb_particles, dt, Tsave = get_info(raw_data_file)

    x_o = Observable(x[1])
    y_o = Observable(y[1])


    
    fig = Figure(y[1])

end