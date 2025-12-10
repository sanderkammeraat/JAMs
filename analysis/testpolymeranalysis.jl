using JLD2
using CairoMakie

CairoMakie.activate!(type = "png")


function load_file(file_location)

    file = jldopen(file_location, "r",iotype=IOStream)

end

#for windows
<<<<<<< HEAD
test_file = load_file(raw"E:\martin\sim_data\p_0.01\raw_data.h5")
=======
#test_file = load_file(raw"E:\martin\sim_data\p_0.1\raw_data.h5")
>>>>>>> 8703c0e46f5957dbebc8bf72632522abe5702e07

#for linux
test_file = load_file("/run/media/martin/HENKESGRFAT/martin/sim_data/p_0.1/raw_data.h5")


function get_info(file)
    
    numb_frames = length(file["frames"])
    numb_particles = length(file["frames"]["1"]["xuw"])


    v_x = zeros(numb_frames, numb_particles)
    v_y = zeros(numb_frames, numb_particles)

    x = zeros(numb_frames, numb_particles)
    y = zeros(numb_frames, numb_particles)

    for i in 1:numb_frames

        x[i,:] = file["frames"]["$i"]["xuw"]
        y[i,:] = file["frames"]["$i"]["yuw"]
        v_x[i,:] = file["frames"]["$i"]["vx"]
        v_y[i,:] = file["frames"]["$i"]["vy"]

    end

    return x, y, v_x, v_y, numb_frames, numb_particles, file["integration_info"]["dt"], file["integration_info"]["Tsave"]
end


function average_velocity(v_x, v_y, numb_frames, numb_particles, sliding_window)
  
    v = zeros(numb_frames-sliding_window+1)

    for j in 1:sliding_window
        v[1] += sum(sqrt.(v_x[j].^2 + v_y[j].^2))/numb_particles
    end

    for i in 2:numb_frames-sliding_window+1
        v[i] += v[i-1] + sum(sqrt.(v_x[i+sliding_window-1].^2 + v_y[i+sliding_window-1].^2))/numb_particles - sum(sqrt.(v_x[i-1].^2 + v_y[i-1].^2))/numb_particles
    end

    return v/sliding_window
end


function plot_avg_velocity!(velocities, numb_frames, numb_particles, sliding_window)

    f = Figure()
    ax = Axis(f[1, 1])

    time = 1:numb_frames-sliding_window+1

    for velocity in velocities
        lines!(time, average_velocity(velocity[1], velocity[2], numb_frames, numb_particles, sliding_window))
    end
    f
end


function MSD(x, y, numb_frames, numb_particles)

    MSD = zeros(numb_frames)

    for Δt in 1:numb_frames
        MSD[Δt] = sum((x[Δt+1:numb_frames,:] - x[1:numb_frames-Δt,:]).^2 + (y[Δt+1:numb_frames,:] - y[1:numb_frames-Δt,:]).^2)/numb_particles/(numb_frames-Δt)
    end

    return MSD

end

function plot_MSD!(positions, numb_frames, numb_particles, save)

    f = Figure()
    ax = Axis(f[1, 1], xscale=log10, yscale=log10)

    time = 1:numb_frames

    for position in positions
        lines!(time, MSD(position[1], position[2], numb_frames, numb_particles))
    end
    f
end


function end_to_end_distance(x, y, numb_frames, numb_particles)

    end_to_end_distance = zeros(numb_frames)
    for i in 1:numb_frames
        end_to_end_distance[i]
    end
end


x, y, v_x, v_y, numb_frames, numb_particles, dt, Tsave = get_info(test_file)

plot_MSD!([[x, y]], numb_frames, numb_particles, false)
plot_avg_velocity!([[v_x, v_y]], numb_frames, numb_particles, 100)

MSD_p_0_1 = MSD(x, y, numb_frames, numb_particles)

f = Figure()
ax = Axis(f[1, 1], xscale=log10, yscale=log10)

time = 1:numb_frames

lines!(ax, time, MSD_p_0_1)
f