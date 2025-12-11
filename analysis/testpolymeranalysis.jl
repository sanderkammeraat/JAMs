using JLD2
using CairoMakie

CairoMakie.activate!(type = "png")


function load_file(file_location)

    file = jldopen(file_location, "r",iotype=IOStream)

end

#for windows
#=
p_0_01 = load_file(raw"E:\martin\sim_data\p_0.01\raw_data.h5")
p_0_03 = load_file(raw"E:\martin\sim_data\p_0.03\raw_data.h5")
p_0_05 = load_file(raw"E:\martin\sim_data\p_0.05\raw_data.h5")
p_0_07 = load_file(raw"E:\martin\sim_data\p_0.07\raw_data.h5")
p_0_1 = load_file(raw"E:\martin\sim_data\p_0.1\raw_data.h5")
p_0_2 = load_file(raw"E:\martin\sim_data\p_0.2\raw_data.h5")
p_0_3 = load_file(raw"E:\martin\sim_data\p_0.3\raw_data.h5")
p_0_4 = load_file(raw"E:\martin\sim_data\p_0.4\raw_data.h5")
=#

#for linux
p_0_08 = load_file("/run/media/martin/HENKESGRFAT/martin/sim_data/p_0.08/raw_data.h5")
p_0_09 = load_file("/run/media/martin/HENKESGRFAT/martin/sim_data/p_0.09/raw_data.h5")
p_0_1 = load_file("/run/media/martin/HENKESGRFAT/martin/sim_data/p_0.1/raw_data.h5")
p_0_11 = load_file("/run/media/martin/HENKESGRFAT/martin/sim_data/p_0.11/raw_data.h5")
p_0_15 = load_file("/run/media/martin/HENKESGRFAT/martin/sim_data/p_0.15/raw_data.h5")


function get_info(file)
    
    numb_frames = length(file["frames"])
    numb_particles = length(file["frames"]["1"]["xuw"])
    pol_id = file["frames/1/pol_id"]
    id_in_pol = file["frames/1/id_in_pol"]


    v_x = zeros(numb_frames, numb_particles)
    v_y = zeros(numb_frames, numb_particles)

    x = zeros(numb_frames, numb_particles)
    y = zeros(numb_frames, numb_particles)

    p_x = zeros(numb_frames, numb_particles)
    p_y = zeros(numb_frames, numb_particles)

    for i in 1:numb_frames

        x[i,:] = file["frames"]["$i"]["xuw"]
        y[i,:] = file["frames"]["$i"]["yuw"]
        v_x[i,:] = file["frames"]["$i"]["vx"]
        v_y[i,:] = file["frames"]["$i"]["vy"]
        p_x[i,:] = file["frames"]["$i"]["px"]
        p_y[i,:] = file["frames"]["$i"]["py"]

    end

    return [x, y, v_x, v_y], p_x, p_y, pol_id, id_in_pol, numb_frames, numb_particles, file["integration_info"]["dt"], file["integration_info"]["Tsave"]
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
    display(f)
end


function MSD(x, y, numb_frames, numb_particles)

    MSD = zeros(convert(Int64, floor((numb_frames-1)/4)))

    for Δt in 1:4:numb_frames-17
        MSD[convert(Int64, floor((Δt-1)/4))+1] += sum((x[Δt+1:numb_frames,:] - x[1:numb_frames-Δt,:]).^2 + (y[Δt+1:numb_frames,:] - y[1:numb_frames-Δt,:]).^2)/numb_particles/(numb_frames-Δt)
    end

    return MSD

end

function plot_MSD!(data, end_time)

    f = Figure()

    ax = Axis(f[1, 1], xscale=log10, yscale=log10)

    lines!(1:end_time/(length(data)-1):end_time, data)
    display(f)
end

function save_MSD(data, folder_name, path)

end


function end_to_end_distance(x, y, pol_id, id_in_pol, numb_frames, numb_particles)

    length_polymer = max(id_in_pol)
    numb_polymers = numb_particles/
    end_to_end_distance = zeros(numb_frames)
    for i in 1:numb_polymers
        end_to_end_distance[i] = 2
    end
end


function radius_of_gyration(x, y, pol_id, id_in_pol, numb_frames, numb_particles)

    length_polymer = max(id_in_pol)
    numb_polymers = numb_particles/length_polymer
    R_2 = zeros(numb_frames)

    center_x = zeros(numb_particles)
    center_y = zeros(numb_particles)

    for i in 1:numb_frames

        for j in 1:numb_polymers
            center_x[(j-1)*length_polymer+1,j*length_polymer] .+= sum(x[(j-1)*length_polymer+1:j*length_polymer])/length_polymer
            center_y[(j-1)*length_polymer+1,j*length_polymer] .+= sum(x[(j-1)*length_polymer+1:j*length_polymer])/length_polymer
        end

        R_2[i] += sum((x - center_x)^2 + (y - center_y)^2)/numb_particles
    end
    return R_2
end

function plot_radius_of_gyration!(multiple, datas, p, end_time)

    if !multiple

        f = Figure()

        ax = Axis(f[1, 1], xscale=log10, yscale=log10)

        lines!(1:end_time/(length(datas)-1):end_time, datas)
        display(f)

    
    else
        f = Figure()
        ax = Axis(f[1, 1])

        i=1
        for data in datas
            activity = p[i]
            lines!(1:end_time/(length(datas)-1):end_time, data, label="p=$activity")
            i+=1
        end
    end

    Legend(f[1,2], ax, "activity")
    display(f)

end


info_0_08, p_x, p_y, pol_id, id_in_pol, numb_frames, numb_particles, dt, Tsave = get_info(p_0_08)
info_0_09, = get_info(p_0_09)
info_0_1, = get_info(p_0_1)
info_0_11, = get_info(p_0_11)
info_0_15, = get_info(p_0_15)


#plot_MSD!([[x, y]], numb_frames, numb_particles, false)
#plot_avg_velocity!([[v_x, v_y]], numb_frames, numb_particles, 100)

MSD_p_0_08 = MSD(info_0_08[1], info_0_08[2], numb_frames, numb_particles)
MSD_p_0_09 = MSD(info_0_09[1], info_0_09[2], numb_frames, numb_particles)
MSD_p_0_1 = MSD(info_0_1[1], info_0_1[2], numb_frames, numb_particles)
MSD_p_0_11 = MSD(info_0_11[1], info_0_11[2], numb_frames, numb_particles)
MSD_p_0_15 = MSD(info_0_15[1], info_0_15[2], numb_frames, numb_particles)


#=
f = Figure()

ax = Axis(f[1, 1], xscale=log10, yscale=log10)

lines!(1:length(MSD_p_0_01), MSD_p_0_01)
f
=#

plot_MSD!(MSD_p_0_08, numb_frames*dt*Tsave)
plot_MSD!(MSD_p_0_09, numb_frames*dt*Tsave)
plot_MSD!(MSD_p_0_1, numb_frames*dt*Tsave)
plot_MSD!(MSD_p_0_11, numb_frames*dt*Tsave)
plot_MSD!(MSD_p_0_15, numb_frames*dt*Tsave)

p = [0.08,0.09,0.1,0.11,0.15]

velocities = [[info_0_01[3], info_0_01[4]],[info_0_03[3], info_0_03[4]],[info_0_05[3], info_0_05[4]],[info_0_07[3], info_0_07[4]],[info_0_1[3], info_0_1[4]],[info_0_2[3], info_0_2[4]],[info_0_3[3], info_0_3[4]],[info_0_4[3], info_0_4[4]]]

plot_avg_velocity!(velocities, p, numb_frames, numb_particles, 800)

R_2_0_08 = radius_of_gyration(info_0_08[1], info_0_08[2], pol_id, id_in_pol, numb_frames, numb_particles)
R_2_0_09 = radius_of_gyration(info_0_09[1], info_0_09[2], pol_id, id_in_pol, numb_frames, numb_particles)
R_2_0_1 = radius_of_gyration(info_0_1[1], info_0_1[2], pol_id, id_in_pol, numb_frames, numb_particles)
R_2_0_11 = radius_of_gyration(info_0_11[1], info_0_11[2], pol_id, id_in_pol, numb_frames, numb_particles)
R_2_0_15 = radius_of_gyration(info_0_15[1], info_0_15[2], pol_id, id_in_pol, numb_frames, numb_particles)

plot_radius_of_gyration!(true, [R_2_0_08,R_2_0_09,R_2_0_1,R_2_0_11,R_2_0_15], p, end_time)