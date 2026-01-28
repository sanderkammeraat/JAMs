using CairoMakie
using JLD2

include("loaddata.jl")

CairoMakie.activate!(type = "pdf")

function plot_data!(data, time, scale, savebool, path, folder_name, file_name, parameters, title, xlabel, ylabel)

    f = Figure()
    ax = Axis(f[1, 1], xscale=scale, yscale=scale, title=title, xlabel=xlabel, ylabel=ylabel)

    if scale == log10
        ylims!(ax, 1e-3,maximum(data)*1.5)
    end


    lines!(ax, time, data)
    display(f)

    if savebool

        save_path = path
        parameters_names = parameters[1]
        parameters_values = parameters[2]
        
        for i in length(parameters)
            parameter_name = parameters_names[i]
            parameter_value = parameters_values[i]
            save_path *= "/"*"$parameter_name-$parameter_value"
            if !isdir(joinpath(path, "plots", "$parameter_name-$parameter_value"))
                mkdir(joinpath(path, "plots", "$parameter_name-$parameter_value"))
            end        
        end
        save(joinpath(save_path, file_name), f)
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

mkdir("/run/media/martin/HENKESGRFAT/martin/sim_data/plots")


sim_folder_name = "sim_data"
#for windows
#path_data = joinpath("E:", "martin", sim_folder_name)

#for linux
path_data = joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name)

ksd = [1.]
kbend = [.3]
kstretch = [1.]
fstretch = [.7]
p = [0.04, 0.06, 0.08, 0.1, 0.13, 0.15, 0.2, 0.4]
kactive = [(-1., 0.), (1., 0.), (0., 1.), (0., -1.), (1/sqrt(2), 1/sqrt(2)),(-1/sqrt(2), 1/sqrt(2)),(1/sqrt(2), -1/sqrt(2)),(-1/sqrt(2), -1/sqrt(2))]
Npol = [150]
N = [1500]

dowesave = true

plotMSD = false
plotbasic_MSD = true
plotaverage_velocity = true
plotradius_of_gyration = false
plotend_to_end_distance = true


for ksd_value in ksd
for kbend_value in kbend
for kstretch_value in kstretch
for fstretch_value in fstretch
for p_value in p
for (kpar_value, kperp_value) in kactive
for Npol_value in Npol
for N_value in N

    parameters = "p_$p_value,kpar_$kpar_value,kperp_$kperp_value"  #"ksd_$ksd_value_kbend_$kbend_value_kstretch_$kstretch_value_fstretch_$f_stretch_value_p_$p_value_kperp_$kperp_value_kpar_$kpar_value_Npol_$Npol_value_N_$N_value"

    analysis = load_file(joinpath(path_data, parameters, "analysis.jld2"))

    if plotbasic_MSD
        mkdir(joinpath(path_data, "plots", "basic_MSD"))
        plot_data!(analysis["basic_MSD"], analysis["basic_MSD_time"], log10, dowesave, path_data, "basic_MSD", "$parameters.pdf", [["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value]], "MSD", "log(time)", "log(MSD)")
    end
    if plotMSD
        mkdir(joinpath(path_data, "plots", "MSD"))
        plot_data!(analysis["MSD"], analysis["MSD_time"], log10, dowesave, path_data, "MSD", "$parameters.pdf", [["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value]],"MSD", "log(time)", "log(MSD)")
    end
    if plotaverage_velocity
        mkdir(joinpath(path_data, "plots", "average_velocity"))
        plot_data!(analysis["average_velocity"], analysis["average_velocity_time"], identity, dowesave, path_data, "average_velocity", "$parameters.pdf", [["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value]],"Average velocity over time", "time", "<v>")
    end
    if plotend_to_end_distance
        mkdir(joinpath(path_data, "plots", "end_to_end_distance"))
        plot_data!(analysis["e_to_e_dist"], analysis["e_to_e_dist_time"], identity, dowesave, path_data, "end_to_end_distance", "$parameters.pdf", [["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value]],"Mean end to end distance", "time", "end to end distance")
    end
    if plotradius_of_gyration
        mkdir(joinpath(path_data, "plots", "radius_of_gyration"))
        plot_data!(analysis["radius_of_gyration"], analysis["radius_of_gyration_time"], identity, dowesave, path_data, "radius_of_gyration", "$parameters.pdf", [["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value]],"Mean radius of gyration", "time", "<R^2>")
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