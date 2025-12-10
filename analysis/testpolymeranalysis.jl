using JLD2
using CairoMakie

CairoMakie.activate!(type = "png")


function load_file(file_location)

    file = jldopen(file_location, "r",iotype=IOStream)

end

#for windows
p_0_01 = load_file(raw"E:\martin\sim_data\p_0.01\raw_data.h5")
p_0_03 = load_file(raw"E:\martin\sim_data\p_0.03\raw_data.h5")
p_0_05 = load_file(raw"E:\martin\sim_data\p_0.05\raw_data.h5")
p_0_07 = load_file(raw"E:\martin\sim_data\p_0.07\raw_data.h5")
p_0_1 = load_file(raw"E:\martin\sim_data\p_0.1\raw_data.h5")
p_0_2 = load_file(raw"E:\martin\sim_data\p_0.2\raw_data.h5")
p_0_3 = load_file(raw"E:\martin\sim_data\p_0.3\raw_data.h5")
p_0_4 = load_file(raw"E:\martin\sim_data\p_0.4\raw_data.h5")

#for linux
#test_file = load_file("/run/media/martin/HENKESGRFAT/martin/sim_data/p_0.1/raw_data.h5")


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

    return [x, y, v_x, v_y], numb_frames, numb_particles, file["integration_info"]["dt"], file["integration_info"]["Tsave"]
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


function plot_avg_velocity!(velocities, p, numb_frames, numb_particles, sliding_window)

    f = Figure()
    ax = Axis(f[1, 1])

    time = 1:numb_frames-sliding_window+1
    i=1
    for velocity in velocities
        activity = p[i]
        lines!(time, average_velocity(velocity[1], velocity[2], numb_frames, numb_particles, sliding_window), label="p=$activity")
        i+=1
    end
    Legend(f[1,2], ax, "activity")
    f
end


function MSD(x, y, numb_frames, numb_particles)

    MSD = zeros(convert(Int64, floor((numb_frames-1)/4)))

    for Δt in 1:4:numb_frames-17
        MSD[convert(Int64, floor((Δt-1)/4))+1] += sum((x[Δt+1:numb_frames,:] - x[1:numb_frames-Δt,:]).^2 + (y[Δt+1:numb_frames,:] - y[1:numb_frames-Δt,:]).^2)/numb_particles/(numb_frames-Δt)
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


info_0_01, numb_frames, numb_particles, dt, Tsave = get_info(p_0_01)
info_0_03, = get_info(p_0_03)
info_0_05, = get_info(p_0_05)
info_0_07, = get_info(p_0_07)
info_0_1, = get_info(p_0_1)
info_0_2, = get_info(p_0_2)
info_0_3, = get_info(p_0_3)
info_0_4, = get_info(p_0_4)


#plot_MSD!([[x, y]], numb_frames, numb_particles, false)
#plot_avg_velocity!([[v_x, v_y]], numb_frames, numb_particles, 100)

MSD_p_0_01 = MSD(info_0_01[1], info_0_01[2], numb_frames, numb_particles)
MSD_p_0_03 = MSD(info_0_03[1], info_0_03[2], numb_frames, numb_particles)
MSD_p_0_05 = MSD(info_0_05[1], info_0_05[2], numb_frames, numb_particles)
MSD_p_0_07 = MSD(info_0_07[1], info_0_07[2], numb_frames, numb_particles)
MSD_p_0_1 = MSD(info_0_1[1], info_0_1[2], numb_frames, numb_particles)
MSD_p_0_2 = MSD(info_0_2[1], info_0_2[2], numb_frames, numb_particles)
MSD_p_0_3 = MSD(info_0_3[1], info_0_3[2], numb_frames, numb_particles)
MSD_p_0_4 = MSD(info_0_4[1], info_0_4[2], numb_frames, numb_particles)

f = Figure()

ax = Axis(f[1, 1], xscale=log10, yscale=log10)

lines!(1:length(MSD_p_0_01), MSD_p_0_01)
f

p = [0.01,0.03,0.05,0.07,0.1,0.2,0.3,0.4]

velocities = [[info_0_01[3], info_0_01[4]],[info_0_03[3], info_0_03[4]],[info_0_05[3], info_0_05[4]],[info_0_07[3], info_0_07[4]],[info_0_1[3], info_0_1[4]],[info_0_2[3], info_0_2[4]],[info_0_3[3], info_0_3[4]],[info_0_4[3], info_0_4[4]]]

plot_avg_velocity!(velocities, p, numb_frames, numb_particles, 800)