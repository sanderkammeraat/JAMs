using CairoMakie
using JLD2

include("loaddata.jl")

CairoMakie.activate!(type = "pdf")

function plot_data!(data, time, scale, savebool, path, file_name, title, xlabel, ylabel)

    f = Figure()
    ax = Axis(f[1, 1], xscale=scale, yscale=scale, title=title, xlabel=xlabel, ylabel=ylabel)

    if scale == log10
        ylims!(ax, 1e-3,maximum(data)*1.5)
    end


    lines!(ax, time, data)
    display(f)

    if savebool
        save(joinpath(path, "plots", file_name), f)
    end
end


function plot_radius_of_gyration!(multiple, datas, p, t_stop)

    if !multiple

        f = Figure()

        ax = Axis(f[1, 1], xscale=log10, yscale=log10)
        time = 1:length(datas)
        lines!(time*t_stop/length(data), datas)
        display(f)

    
    else
        f = Figure()
        ax = Axis(f[1, 1])

        i=1
        for data in datas
            activity = p[i]
            time = 1:length(data)
            lines!(time*t_stop/length(data), data, label="p=$activity")
            i+=1
        end
    end

    Legend(f[1,2], ax, "activity")
    display(f)

end


sim_folder_name = "sim_data"
#for windows
#path_data = joinpath("E:", "martin", sim_folder_name)

#for linux
path_data = joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name)

ksd = [1.]
kbend = [.3]
kstretch = [1.]
fstretch = [.7]
p = [0.01, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2]
kperp = [0]
kpar = [-1]
Npol = [150]
N = [1500]

dowesave = true

plotMSD = true
plotaverage_velocity = true
plotradius_of_gyration = false
plotend_to_end_distance = true


for ksd_value in ksd
for kbend_value in kbend
for kstretch_value in kstretch
for fstretch_value in fstretch
for p_value in p
for kperp_value in kperp
for kpar_value in kpar
for Npol_value in Npol
for N_value in N 

    parameters = "p_$p_value" #"ksd_$ksd_value_kbend_$kbend_value_kstretch_$kstretch_value_fstretch_$f_stretch_value_p_$p_value_kperp_$kperp_value_kpar_$kpar_value_Npol_$Npol_value_N_$N_value"

    analysis = load_file(joinpath(path_data, parameters, "analysis.jld2"))

    if plotMSD
        MSD_data = analysis["MSD"][1]
        plot_data!(MSD_data, MSD_time, log10, dowesave, path_data, "MSD/$parameters.png")
    end
    if plotaverage_velocity
        average_velocity_data = analysis["average_velocity"][1]
        plot_data!(average_velocity_data, average_velocity_time, identity, dowesave, path_data, "average_velocity/$parameters.png")
    end
    if plotend_to_end_distance
        end_to_end_distance_data = analysis["end_to_end_distance"][1]
        plot_data!(end_to_end_distance_data, end_to_end_distance_time, identity, dowesave, path_data, "end_to_end_distance/$parameters.png")
    end
    if plotradius_of_gyration
        radius_of_gyration_data = analysis["radius_of_gyration"][1]
        plot_data!(radius_of_gyration_data, radius_of_gyration_time, identity, dowesave, path_data, "radius_of_gyration/$parameters.png")
    end

    close(analysis)
end
end
end
end
end
end
end
end
end