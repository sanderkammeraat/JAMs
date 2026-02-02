using CairoMakie
using JLD2

include("loaddata.jl")


function plot_data!(data, time, scale, displaybool, savebool, path, folder_name, file_name, para_names, para_vals, title, xlabel, ylabel)

    f = Figure()
    ax = Axis(f[1, 1], xscale=scale, yscale=scale, title=title, xlabel=xlabel, ylabel=ylabel)

    if scale == log10
        ylims!(ax, 1e-3,maximum(data)*1.5)
    end

    lines!(ax, time, data)

    if displaybool
        display(f)
    end

    if savebool

        save_path = joinpath(path, folder_name)
        endpoint = length(para_names)
        
        for i in 1:endpoint
            parameter_name = para_names[i]
            parameter_value = para_vals[i]
            save_path = joinpath(save_path, "$parameter_name-$parameter_value")
            if !isdir(save_path)
                mkdir(save_path)
            end        
        end

        save(joinpath(save_path, file_name), f)

    end
end


sim_folder_name = "sim_data"
#for windows
path_data = joinpath("E:", "martin", sim_folder_name)

#for linux
#path_data = joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name)


if !isdir(joinpath(path_data, "plots"))
    mkdir(joinpath(path_data, "plots"))
end
if !isdir(joinpath(path_data, "plots", "basic_MSD"))
    mkdir(joinpath(path_data, "plots", "basic_MSD"))
end
if !isdir(joinpath(path_data, "plots", "MSD"))
    mkdir(joinpath(path_data, "plots", "MSD"))
end
if !isdir(joinpath(path_data, "plots", "average_velocity"))
    mkdir(joinpath(path_data, "plots", "average_velocity"))
end
if !isdir(joinpath(path_data, "plots", "end_to_end_distance"))
    mkdir(joinpath(path_data, "plots", "end_to_end_distance"))
end
if !isdir(joinpath(path_data, "plots", "radius_of_gyration"))
    mkdir(joinpath(path_data, "plots", "radius_of_gyration"))
end


ksd = [1.]
kbend = [.3]
kstretch = [1.]
fstretch = [.7]
p = [0.04, 0.06, 0.08, 0.1, 0.13, 0.15, 0.2, 0.4]
kactive = [(-1., 0.), (1., 0.), (0., 1.), (0., -1.), (1/sqrt(2), 1/sqrt(2)),(-1/sqrt(2), 1/sqrt(2)),(1/sqrt(2), -1/sqrt(2)),(-1/sqrt(2), -1/sqrt(2))]
Npol = [150]
N = [1500]

dowesave = true
doweplot = false

plotMSD = true
plotbasic_MSD = false
plotaverage_velocity = false
plotradius_of_gyration = false
plotend_to_end_distance = false


for ksd_value in ksd
for kbend_value in kbend
for kstretch_value in kstretch
for fstretch_value in fstretch
for p_value in p
for (kpar_value, kperp_value) in kactive
for Npol_value in Npol
for N_value in N

    param = "p_$p_value,kpar_$kpar_value,kperp_$kperp_value"  #"ksd_$ksd_value_kbend_$kbend_value_kstretch_$kstretch_value_fstretch_$f_stretch_value_p_$p_value_kperp_$kperp_value_kpar_$kpar_value_Npol_$Npol_value_N_$N_value"

    analysis = load_file(joinpath(path_data, param, "analysis.jld2"))

    if plotbasic_MSD
        plot_data!(analysis["basic_MSD"], analysis["basic_MSD_time"], log10, doweplot, dowesave, joinpath(path_data, "plots"), "basic_MSD", "$param.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value], "MSD", "log(time)", "log(MSD)")
    end
    if plotMSD
        plot_data!(analysis["MSD"], analysis["MSD_time"], log10, doweplot, dowesave, joinpath(path_data, "plots"), "MSD", "$param.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"MSD", "log(time)", "log(MSD)")
    end
    if plotaverage_velocity
        plot_data!(analysis["average_velocity"], analysis["average_velocity_time"], identity, doweplot, dowesave, joinpath(path_data, "plots"), "average_velocity", "$param.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Average velocity over time", "time", "<v>")
    end
    if plotend_to_end_distance
        plot_data!(analysis["e_to_e_dist"], analysis["e_to_e_dist_time"], identity, doweplot, dowesave, joinpath(path_data, "plots"), "end_to_end_distance", "$param.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Mean end to end distance", "time", "end to end distance")
    end
    if plotradius_of_gyration
        plot_data!(analysis["radius_of_gyration"], analysis["radius_of_gyration_time"], identity, doweplot, dowesave, joinpath(path_data, "plots"), "radius_of_gyration", "$param.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Mean radius of gyration", "time", "<R^2>")
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