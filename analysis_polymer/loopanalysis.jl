include("analysisfunctions.jl")
include("plotanalysis.jl")
include("replotpolymer.jl")
include("loaddata.jl")


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
kactive = [(-1., 0.), (1., 0.)#=, (0., 1.), (0., -1.), (1/sqrt(2), 1/sqrt(2)),(-1/sqrt(2), 1/sqrt(2)),(1/sqrt(2), -1/sqrt(2)),(-1/sqrt(2), -1/sqrt(2))=#]
Npol = [750]
N = [2250]


MSDbool = false
polymer_MSDbool = false
basic_MSDbool = false
average_velocitybool = true
radius_of_gyrationbool = true
end_to_end_distancebool = true


save_plotbool = true
replotbool = false

if save_plotbool
    create_directory(path_data, "plots")
end


for ksd_value in ksd
for kbend_value in kbend
for kstretch_value in kstretch
for fstretch_value in fstretch
for (i, p_value) in enumerate(p)
for (j, (kpar_value, kperp_value)) in enumerate(kactive)
for Npol_value in Npol
for N_value in N

    parameters = "p_$p_value,kpar_$kpar_value,kperp_$kperp_value" #"ksd_$ksd_value_kbend_$kbend_value_kstretch_$kstretch_value_fstretch_$f_stretch_value_p_$p_value_kperp_$kperp_value_kpar_$kpar_value_Npol_$Npol_value_N_$N_value"

    max = length(p)*length(kactive)
    progress = convert(Int64, floor(100*((i-1)*length(kactive) + j-1)/max))
    println("$parameters\n\n$progress% done\n\n")


    dataset = get_data(joinpath(path_data, "$parameters", "raw_data.h5"))

    if isfile(joinpath(path_data, "$parameters", "analysis.jld2"))

        println("\narrived in loop where analysis already exists\n")

        data = jldopen(joinpath(path_data, "$parameters", "analysis.jld2"), "w")


        if MSDbool

            MSD_data = MSD(dataset.x, dataset.y, dataset.numb_frames, dataset.N, dataset.t_stop)
            data["MSD"] = MSD_data["MSD"]
            data["MSD_time"] = MSD_data["MSD_time"]

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "MSD")
                plot_data!(MSD_data["MSD"], MSD_data["MSD_time"], log10, false, true, joinpath(path_data, "plots"), "MSD", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"MSD", "log(time)", "log(MSD)")
            end
        end


        if polymer_MSDbool

            polymer_MSD_data = polymer_MSD(dataset.x, dataset.y, dataset.pol_id, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop)
            data["polymer_MSD"] = MSD_data["polymer_MSD"]
            data["polymer_MSD_time"] = MSD_data["polymer_MSD_time"]

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "polymer_MSD")
                plot_data!(polymer_MSD_data["polymer_MSD"], polymer_MSD_data["polymer_MSD_time"], log10, false, true, joinpath(path_data, "plots"), "polymer_MSD", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"MSD of the center of mass of the polymers", "log(time)", "log(MSD)")
            end
        end


        if basic_MSDbool

            basic_MSD_data = basic_MSD(dataset.x, dataset.y, dataset.numb_frames, dataset.N, dataset.t_stop)
            data["basic_MSD"] = basic_MSD_data["basic_MSD"]
            data["basic_MSD_time"] = basic_MSD_data["basic_MSD_time"]

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "basic_MSD")
                plot_data!(basic_MSD_data["basic_MSD"], basic_MSD_data["basic_MSD_time"], log10, false, true, joinpath(path_data, "plots"), "MSD", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"MSD", "log(time)", "log(MSD)")
            end
        end


        if average_velocitybool

            average_velocity_data = average_velocity(dataset.vx, dataset.vy, dataset.numb_frames, dataset.N, dataset.t_stop)
            data["average_velocity"] = average_velocity_data["average_velocity"]
            data["average_velocity_time"] = average_velocity_data["average_velocity_time"]

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "average_velocity")
                plot_data!(average_velocity_data["average_velocity"], average_velocity_data["average_velocity_time"], log10, false, true, joinpath(path_data, "plots"), "average_velocity", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Average velocity over time", "time", "<v>")
            end
        end


        if radius_of_gyrationbool
            radius_of_gyration_data = radius_of_gyration(dataset.x, dataset.y, dataset.pol_id, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop)
            data["R_2"] = radius_of_gyration_data["R_2"]
            data["R_2_time"] = radius_of_gyration_data["R_2_time"]

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "R_2")
                plot_data!(radius_of_gyration_data["R_2"], radius_of_gyration_data["R_2_time"], log10, false, true, joinpath(path_data, "plots"), "radius_of_gyration", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Mean radius of gyration", "time", "<R^2>")
            end
        end


        if end_to_end_distancebool

            end_to_end_distance_data = end_to_end_distance(dataset.x, dataset.y, dataset.pol_id, dataset.id_in_pol, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop)
            data["e_to_e_dist"] = end_to_end_distance_data["e_to_e_dist"]
            data["e_to_e_dist_time"] = end_to_end_distance_data["e_to_e_dist_time"]

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "e_to_e_dist")
                plot_data!(end_to_end_distance_data["e_to_e_dist"], end_to_end_distance_data["e_to_e_dist_time"], identity, false, true, joinpath(path_data, "plots"), "end_to_end_distance", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Mean end to end distance", "time", "end to end distance")
            end
        end

        close(data)

    else

        data = Dict{String, Vector{Float64}}()
        println("\nanalysis does not exist")

        if MSDbool

            MSD_data = MSD(dataset.x, dataset.y, dataset.numb_frames, dataset.N, dataset.t_stop)
            merge!(data, MSD_data)

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "MSD")
                plot_data!(MSD_data["MSD"], MSD_data["MSD_time"], log10, false, true, joinpath(path_data, "plots"), "MSD", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"MSD", "log(time)", "log(MSD)")
            end
        end


        if polymer_MSDbool

            polymer_MSD_data = polymer_MSD(dataset.x, dataset.y, dataset.pol_id, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop)
            data["polymer_MSD"] = MSD_data["polymer_MSD"]
            data["polymer_MSD_time"] = MSD_data["polymer_MSD_time"]

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "polymer_MSD")
                plot_data!(polymer_MSD_data["polymer_MSD"], polymer_MSD_data["polymer_MSD_time"], log10, false, true, joinpath(path_data, "plots"), "polymer_MSD", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"MSD of the center of mass of the polymers", "log(time)", "log(MSD)")
            end
        end


        if basic_MSDbool

            basic_MSD_data = basic_MSD(dataset.x, dataset.y, dataset.numb_frames, dataset.N, dataset.t_stop)
            merge!(data, basic_MSD_data)

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "basic_MSD")
                plot_data!(basic_MSD_data["basic_MSD"], basic_MSD_data["basic_MSD_time"], log10, false, true, joinpath(path_data, "plots"), "MSD", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"MSD", "log(time)", "log(MSD)")
            end
        end


        if average_velocitybool

            println("\narrived in avg vel from non preexisting analysis")

            average_velocity_data = average_velocity(dataset.vx, dataset.vy, dataset.numb_frames, dataset.N, dataset.t_stop)
            merge!(data, average_velocity_data)

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "average_velocity")
                plot_data!(average_velocity_data["average_velocity"], average_velocity_data["average_velocity_time"], log10, false, true, joinpath(path_data, "plots"), "average_velocity", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Average velocity over time", "time", "<v>")
            end
        end


        if radius_of_gyrationbool

            println("\narrived in radius gyr from non preexisting analysis")

            radius_of_gyration_data = radius_of_gyration(dataset.x, dataset.y, dataset.pol_id, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop)
            merge!(data, radius_of_gyration_data)

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "R_2")
                plot_data!(radius_of_gyration_data["R_2"], radius_of_gyration_data["R_2_time"], log10, false, true, joinpath(path_data, "plots"), "radius_of_gyration", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Mean radius of gyration", "time", "<R^2>")
            end
        end


        if end_to_end_distancebool

            println("\narrived in end_to_end_distance from non preexisting analysis")

            end_to_end_distance_data = end_to_end_distance(dataset.x, dataset.y, dataset.pol_id, dataset.id_in_pol, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop)
            merge!(data, end_to_end_distance_data)

            if save_plotbool
                create_directory(joinpath(path_data, "plots"), "e_to_e_dist")
                plot_data!(end_to_end_distance_data["e_to_e_dist"], end_to_end_distance_data["e_to_e_dist_time"], identity, false, true, joinpath(path_data, "plots"), "end_to_end_distance", "$parameters.pdf", ["p", "kpar", "kperp"], [p_value, kpar_value, kperp_value],"Mean end to end distance", "time", "end to end distance")
            end
        end

        save_data(data, joinpath(path_data, "$parameters"))
        
    end

end
end
end
end
end
end
end
end


if replotbool
    for p_value in p
    for (kpar_value, kperp_value) in kactive

        param = "p_$p_value,kpar_$kpar_value,kperp_$kperp_value"

        if !isfile(joinpath(path_data, "$param", "movie", "sim_movie.mp4"))

            raw_data_file = load_file(joinpath(path_data, "$param", "raw_data.h5"))

            make_movie(raw_data_file,joinpath(path_data, "$param", "movie"))

            close(raw_data_file)
            
            GLMakie.closeall()
            
        end
        
    end
    end
end