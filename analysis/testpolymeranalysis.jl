using JLD2
using CairoMakie


function load_file(file_location)

    file = jldopen(file_location, "r",iotype=IOStream)

end

#for windows
test_file = load_file(raw"E:\martin\sim_data\p_0.1\raw_data.h5")

#for linux
#test_file = load_file("/run/media/martin/HENKESGRFAT/martin/sim_data/p_0.1/raw_data.h5")


function average_velocity(file, sliding_window)


    numb_frames = length(file["frames"])
    numb_particles = length(file["frames"]["1"]["x"])


    v_x = zeros(numb_frames, numb_particles)
    v_y = zeros(numb_frames, numb_particles)

    for i in 1:numb_frames

        v_x[i,:] = file["frames"]["$i"]["vx"]
        v_y[i,:] = file["frames"]["$i"]["vy"]

    end

    v = zeros(numb_frames-sliding_window+1)

    for j in 1:sliding_window
        v[1] += sum(sqrt.(v_x[j].^2 + v_y[j].^2))/numb_particles
    end

    for i in 2:numb_frames-sliding_window+1
        v[i] += v[i-1] + sum(sqrt.(v_x[i+sliding_window-1].^2 + v_y[i+sliding_window-1].^2))/numb_particles - sum(sqrt.(v_x[i-1].^2 + v_y[i-1].^2))/numb_particles
    end

    return v/sliding_window
end


function plot_avg_velocity!(files, sliding_window)

    f = Figure()
    ax = Axis(f[1, 1])

    delta_t = files[1]["integration_info"]["dt"]
    t_stop = files[1]["integration_info"]["t_stop"]
    Tsave = files[1]["integration_info"]["Tsave"]
    numb_frames = length(files[1]["frames"])

    time = 1:numb_frames-sliding_window+1

    for file in files
        lines!(time, average_velocity(file, sliding_window))
    end
    f
end


function MSD(file, sliding_window)

    numb_frames = length(file["frames"])
    numb_particles = length(file["frames"]["1"]["x"])


    x = zeros(numb_frames, numb_particles)
    y = zeros(numb_frames, numb_particles)

    for i in 1:numb_frames

        x[i,:] = file["frames"]["$i"]["xuw"]
        y[i,:] = file["frames"]["$i"]["yuw"]

    end

    MSD = zeros(numb_frames)

    for i in 1:numb_frames-sliding_window+1
        for j in i:i+sliding_window-1
            MSD[i] += sum((x[j]-x[i]).^2 .+ (y[j]-y[i]).^2)/numb_particles
        end
    end

    return MSD/sliding_window

end

function plot_MSD!(files, sliding_window)

    f = Figure()
    ax = Axis(f[1, 1], xscale=log10, yscale=log10)
    ylims!(ax, 1e-6, 1e-2)
    delta_t = files[1]["integration_info"]["dt"]
    t_stop = files[1]["integration_info"]["t_stop"]
    Tsave = files[1]["integration_info"]["Tsave"]
    numb_frames = length(files[1]["frames"])

    time = 1:numb_frames

    for file in files
        lines!(time, MSD(file, sliding_window))
    end
    f
end


plot_MSD!([test_file], 100)
plot_avg_velocity!([test_file], 100)