using CairoMakie
using JLD2

CairoMakie.activate!(type = "png")



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


function plot_MSD!(data, end_time)

    f = Figure()
    ax = Axis(f[1, 1], xscale=log10, yscale=log10)
    time = 1:length(data)
    ylims!(ax, 1e-3,maximum(data)*1.5)


    lines!(ax, time*end_time/length(data), data)
    display(f)
end


function plot_radius_of_gyration!(multiple, datas, p, end_time)

    if !multiple

        f = Figure()

        ax = Axis(f[1, 1], xscale=log10, yscale=log10)
        time = 1:length(datas)
        lines!(time*end_time/length(data), datas)
        display(f)

    
    else
        f = Figure()
        ax = Axis(f[1, 1])

        i=1
        for data in datas
            activity = p[i]
            time = 1:length(data)
            lines!(time*end_time/length(data), data, label="p=$activity")
            i+=1
        end
    end

    Legend(f[1,2], ax, "activity")
    display(f)

end


#=
info, p_x, p_y, pol_id, id_in_pol, numb_frames, numb_particles, dt, Tsave = get_info(file)

MSD_file = MSD(info[1], info[2], numb_frames, numb_particles)
average_velocity_file = average_velocity(info[3], info[4], numb_frames, numb_particles, 2)
radius_of_gyration_file = radius_of_gyration(info[1], info[2], pol_id, id_in_pol, numb_frames, numb_particles)

plot_MSD!(MSD_file, numb_frames*dt*Tsave)
save_data(MSD_file, "MSD", raw"C:\Users\gabri\Documents\Travail-Etude\Master's Theoretical Physics Leiden\Research Project\Data\sim_data\p_0.13")

saved_data = jldopen("C:\\Users\\gabri\\Documents\\Travail-Etude\\Master's Theoretical Physics Leiden\\Research Project\\JAMs\\MSD.jld2", "r")

plot_MSD!(saved_data["data"], numb_frames*dt*Tsave)

=#


info_0_01, p_x, p_y, pol_id, id_in_pol, numb_frames, numb_particles, dt, Tsave = get_info(p_0_01)
info_0_04, = get_info(p_0_04)
info_0_06, = get_info(p_0_06)
info_0_08, = get_info(p_0_08)
info_0_1, = get_info(p_0_1)
info_0_15, = get_info(p_0_15)
info_0_2, = get_info(p_0_2)


close(p_0_01)
close(p_0_04)
close(p_0_06)
close(p_0_08)
close(p_0_1)
close(p_0_15)
close(p_0_2)
close(p_0_3)

MSD_p_0_01 = MSD(info_0_01[1], info_0_01[2], numb_frames, numb_particles)
MSD_p_0_04 = MSD(info_0_04[1], info_0_04[2], numb_frames, numb_particles)
MSD_p_0_06 = MSD(info_0_06[1], info_0_06[2], numb_frames, numb_particles)
MSD_p_0_08 = MSD(info_0_08[1], info_0_08[2], numb_frames, numb_particles)
MSD_p_0_1 = MSD(info_0_1[1], info_0_1[2], numb_frames, numb_particles)
MSD_p_0_15 = MSD(info_0_15[1], info_0_15[2], numb_frames, numb_particles)
MSD_p_0_2 = MSD(info_0_2[1], info_0_2[2], numb_frames, numb_particles)

#=
plot_MSD!(MSD_p_0_08, numb_frames*dt*Tsave)
plot_MSD!(MSD_p_0_09, numb_frames*dt*Tsave)
plot_MSD!(MSD_p_0_1, numb_frames*dt*Tsave)
plot_MSD!(MSD_p_0_11, numb_frames*dt*Tsave)
plot_MSD!(MSD_p_0_15, numb_frames*dt*Tsave)
=#



velocities = [[info_0_01[3], info_0_01[4]],[info_0_04[3], info_0_04[4]],[info_0_06[3], info_0_06[4]],[info_0_08[3], info_0_08[4]],[info_0_1[3], info_0_1[4]],[info_0_15[3], info_0_15[4]],[info_0_2[3], info_0_2[4]]]


save_data([MSD_p_0_01, velocities[1]], ["MSD", "average_velocities"], raw"E:\martin\sim_data\p_0.01")
save_data([MSD_p_0_04, velocities[2]], ["MSD", "average_velocities"], raw"E:\martin\sim_data\p_0.04")
save_data([MSD_p_0_06, velocities[3]], ["MSD", "average_velocities"], raw"E:\martin\sim_data\p_0.06")
save_data([MSD_p_0_08, velocities[4]], ["MSD", "average_velocities"], raw"E:\martin\sim_data\p_0.08")
save_data([MSD_p_0_1, velocities[5]], ["MSD", "average_velocities"], raw"E:\martin\sim_data\p_0.1")
save_data([MSD_p_0_15, velocities[6]], ["MSD", "average_velocities"], raw"E:\martin\sim_data\p_0.15")
save_data([MSD_p_0_2, velocities[7]], ["MSD", "average_velocities"], raw"E:\martin\sim_data\p_0.2")
save_data([MSD_p_0_3, velocities[8]], ["MSD", "average_velocities"], raw"E:\martin\sim_data\p_0.3")

#plot_avg_velocity!(velocities, p, numb_frames, numb_particles, 1000)

#=
R_2_0_08 = radius_of_gyration(info_0_08[1], info_0_08[2], pol_id, id_in_pol, numb_frames, numb_particles)
R_2_0_09 = radius_of_gyration(info_0_09[1], info_0_09[2], pol_id, id_in_pol, numb_frames, numb_particles)
R_2_0_1 = radius_of_gyration(info_0_1[1], info_0_1[2], pol_id, id_in_pol, numb_frames, numb_particles)
R_2_0_11 = radius_of_gyration(info_0_11[1], info_0_11[2], pol_id, id_in_pol, numb_frames, numb_particles)
R_2_0_15 = radius_of_gyration(info_0_15[1], info_0_15[2], pol_id, id_in_pol, numb_frames, numb_particles)

convert(Array{Float64},)

plot_radius_of_gyration!(true, [vec(R_2_0_08),vec(R_2_0_09),vec(R_2_0_1),vec(R_2_0_11),vec(R_2_0_15)], p, numb_frames*dt*Tsave)

=#