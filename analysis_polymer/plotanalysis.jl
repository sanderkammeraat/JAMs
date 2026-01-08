using CairoMakie
using JLD2

include("loopanalysis.jl")

CairoMakie.activate!(type = "png")

function plot_data!(data, scale, end_time, savebool)

    f = Figure()
    ax = Axis(f[1, 1], xscale=scale, yscale=scale)
    time = 1:length(data)
    #ylims!(ax, 1e-3,maximum(data)*1.5)
    display(typeof(data))

    lines!(ax, time*end_time/length(data), data)
    display(f)
    if savebool
        save(raw"E:\martin\sim_data\plots\average_velocity\p_0_01.png", f)
    end
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


file = load_file(raw"E:\martin\sim_data\p_0.01\analysis.jld2")
average_velocity_data = file["average_velocity"][1]
close(file)

plot_data!(average_velocity_data, identity, 51, true)


dowesave = true

plotMSD = true
plotaverage_velocity = true
plotradius_of_gyration = true
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

    #for windows
    path_data = joinpath("E:", "martin", sim_folder_name, parameters)

    #for linux
    #path_data = joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name, parameters)

    analysis = load_file(joinpath(path_data, parameters, "analysis.jld2"))

    if plotMSD & MSDbool
        MSD_data = analysis["MSD"][1]
        plot_data!(MSD_data, log10, 51, dowesave)
    end
    if plotaverage_velocity & average_velocitybool
    end
    if plotend_to_end_distance & end_to_end_distancebool
    end
    if plotradius_of_gyration & radius_of_gyrationbool
    end
end
end
end
end
end
end
end
end
end