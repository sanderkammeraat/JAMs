using CairoMakie
using JLD2

if false
    include("loaddata.jl")

    CairoMakie.activate!()

    sim_folder_name = "sim_data"

    #for windows
    #path_data = joinpath("E:", "martin", sim_folder_name)

    #for linux
    #path_data = joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name)

    path_data = joinpath(raw"C:\Users\gabri\Documents\Travail-Etude\Master's Theoretical Physics Leiden\Research Project\Data", "sim_data")


end


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


function create_directory(path_to_folder, folder_name)

    if !isdir(joinpath(path_to_folder, folder_name))
        mkdir(joinpath(path_to_folder, folder_name))
    end

end

if false

    ksd = [1.]
    kbend = [.3]
    kstretch = [1.]
    fstretch = [.7]
    p = [0.04, 0.06, 0.08, 0.1, 0.13, 0.15, 0.2, 0.4]
    kactive = [(-1., 0.), (1., 0.)#=, (0., 1.), (0., -1.), (1/sqrt(2), 1/sqrt(2)),(-1/sqrt(2), 1/sqrt(2)),(1/sqrt(2), -1/sqrt(2)),(-1/sqrt(2), -1/sqrt(2))=#]
    Npol = [150]
    N = [1500]

    dowesave = true
    doweplot = false

    plotMSD = false
    plotpolymer_MSD = false
    plotbasic_MSD = false
    plotaverage_velocity = true
    plotradius_of_gyration = true
    plotend_to_end_distance = true


    create_directory(path_data, "plots")


    for ksd_value in ksd
    for kbend_value in kbend
    for kstretch_value in kstretch
    for fstretch_value in fstretch
    for (i, p_value) in enumerate(p)
    for (j, (kpar_value, kperp_value)) in enumerate(kactive)
    for Npol_value in Npol
    for N_value in N

        parameters = "p_$p_value,kpar_$kpar_value,kperp_$kperp_value"  #"ksd_$ksd_value_kbend_$kbend_value_kstretch_$kstretch_value_fstretch_$f_stretch_value_p_$p_value_kperp_$kperp_value_kpar_$kpar_value_Npol_$Npol_value_N_$N_value"

        
        max = length(p)*length(kactive)
        progress = convert(Int64, floor(100*((i-1)*length(kactive) + j-1)/max))
        println("$parameters\n\n$progress% done\n\n")

        
        analysis = load_file(joinpath(path_data, parameters, "analysis.jld2"))

        if plotbasic_MSD
            create_directory(joinpath(path_data, "plots"), "basic_MSD")
            plot_data!(analysis["basic_MSD"], analysis["basic_MSD_time"], log10, doweplot, dowesave, joinpath(path_data, "plots"), "basic_MSD", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value], "MSD", "log(time)", "log(MSD)")
        end


        if plotMSD
            create_directory(joinpath(path_data, "plots"), "MSD")
            plot_data!(analysis["MSD"], analysis["MSD_time"], log10, doweplot, dowesave, joinpath(path_data, "plots"), "MSD", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"MSD", "log(time)", "log(MSD)")
        end


        if plotpolymer_MSD
            create_directory(joinpath(path_data, "plots"), "polymer_MSD")
            plot_data!(analysis["polymer_MSD"], analysis["polymer_MSD_time"], log10, doweplot, dowesave, joinpath(path_data, "plots"), "polymer_MSD", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"MSD of the center of mass of the polymers", "log(time)", "log(MSD)")
        end


        if plotaverage_velocity
            create_directory(joinpath(path_data, "plots"), "average_velocity")
            plot_data!(analysis["average_velocity"], analysis["average_velocity_time"], identity, doweplot, dowesave, joinpath(path_data, "plots"), "average_velocity", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Average velocity over time", "time", "<v>")
        end


        if plotend_to_end_distance
            create_directory(joinpath(path_data, "plots"), "e_to_e_dist")
            plot_data!(analysis["e_to_e_dist"], analysis["e_to_e_dist_time"], identity, doweplot, dowesave, joinpath(path_data, "plots"), "end_to_end_distance", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Mean end to end distance", "time", "end to end distance")
        end


        if plotradius_of_gyration
            create_directory(joinpath(path_data, "plots"), "R_2")
            plot_data!(analysis["R_2"], analysis["R_2_time"], identity, doweplot, dowesave, joinpath(path_data, "plots"), "radius_of_gyration", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Mean radius of gyration", "time", "<R^2>")
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