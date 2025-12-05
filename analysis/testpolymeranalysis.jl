using JLD2
using GLMakie
using CairoMakie


function load_file(file_location)

    file = jldopen(file_location, "r",iotype=IOStream)

end

#=
p_01 = load_file(raw"C:\Users\gabri\Documents\Travail-Etude\Master's Theoretical Physics Leiden\Research Project\Data\test_saving\simdata\p_0.1\raw_data.h5")
p_05 = load_file(raw"C:\Users\gabri\Documents\Travail-Etude\Master's Theoretical Physics Leiden\Research Project\Data\test_saving\simdata\p_0.5\raw_data.h5")
p_10 = load_file(raw"C:\Users\gabri\Documents\Travail-Etude\Master's Theoretical Physics Leiden\Research Project\Data\test_saving\simdata\p_1.0\raw_data.h5")
=#


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

    for i in 1:numb_frames-sliding_window+1
        for j in i:i+sliding_window-1
            v[i] += sum(sqrt.(v_x[j].^2 + v_y[j].^2))/numb_particles
        end
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


function MSD(file)

    numb_frames = length(file["frames"])
    numb_particles = length(file["frames"]["1"]["x"])


    x = zeros(numb_frames, numb_particles)
    y = zeros(numb_frames, numb_particles)

    for i in 1:numb_frames

        x[i,:] = file["frames"]["$i"]["xuw"]
        y[i,:] = file["frames"]["$i"]["yuw"]

    end

    MSD = zeros(numb_frames)

    for i in 1:numb_frames
        
        MSD[i] += sum((x[i]-x[1]).^2 .+ (y[i]-y[1]).^2)/numb_particles
    
    end

    return MSD

end

function plot_MSD!(files)

    f = Figure()
    ax = Axis(f[1, 1])

    delta_t = files[1]["integration_info"]["dt"]
    t_stop = files[1]["integration_info"]["t_stop"]
    Tsave = files[1]["integration_info"]["Tsave"]
    numb_frames = length(files[1]["frames"])

    time = 1:numb_frames

    for file in files
        lines!(time, MSD(file))
    end
    f
end


function isjammed(file, window_percentage, tolerance_percentage, numb_accepted_deviation)
    MSD = MSD(file)
    velocity = average_velocity(file, 1)
    start = length(MSD)-floor(window_percentage/100*length(MSD))
    finish = length(MSD)

    mean_MSD = sum(MSD[start:finish])/(finish-start)
    mean_velocity = sum(velocity[start:finish])/(finish-start)

    count_MSD = 0
    count_veocity = 0

    for i in length(MSD)-floor(window_percentage/100*length(MSD)):length(MSD)
        if abs(MSD[i]-mean_MSD) > tolerance/100*mean_MSD
            count_MSD += 1
        end

        if abs(velocity[i]-mean_velocity) > tolerance/100*mean_velocity
            count_veocity += 1
        end
    end

    if count_MSD > numb_accepted_deviation
        MSD_jammed = true
    else 
        MSD_jammed = false
    end

    if count_velocity > numb_accepted_deviation
        velocity_jammed = true
    else 
        velocity_jammed = false
    end

    return MSD_jammed, velocity_jammed
end

plot_MSD!([p_01, p_05, p_10])
#plot_avg_velocity!([p_01, p_05, p_10], 20)