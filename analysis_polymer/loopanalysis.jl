include("analysisfunctions.jl")
include("loaddata.jl")


sim_folder_name = "sim_data"

#for windows
path_data = joinpath("E:", "martin", sim_folder_name)

#for linux
#path_data = joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name)


ksd = [1.]
kbend = [.3]
kstretch = [1.]
fstretch = [.7]
p = [0., 0.04,#= 0.06, 0.08, =#0.1, #=0.13, 0.15, 0.2, =#0.4]
kpar = [-1]
Npol = [750]
N_in_pol = [2, 10, 20]


MSDbool = false
polymer_MSDbool = false
basic_MSDbool = false
average_velocitybool = false
radius_of_gyrationbool = false
end_to_end_distancebool = true



for ksd_value in ksd
for kbend_value in kbend
for kstretch_value in kstretch
for fstretch_value in fstretch
for (i, p_value) in enumerate(p)
for kpar_value in kpar
for Npol_value in Npol
for (j, N_in_pol_value) in enumerate(N_in_pol)

    parameters = "p_$p_value,N_$N_in_pol_value" #"ksd_$ksd_value_kbend_$kbend_value_kstretch_$kstretch_value_fstretch_$f_stretch_value_p_$p_value_kperp_$kperp_value_kpar_$kpar_value_Npol_$Npol_value_N_$N_value"

    max = length(p)*length(N_in_pol)
    progress = convert(Int64, floor(100*((i-1)*length(N_in_pol) + j-1)/max))
    println("$parameters\n\n$progress% done\n\n")


    dataset = get_data(joinpath(path_data, "$parameters", "raw_data.h5"))

    if isfile(joinpath(path_data, "$parameters", "analysis.jld2"))

        data = jldopen(joinpath(path_data, "$parameters", "analysis.jld2"), "w+")


        if MSDbool

            MSD_data = MSD(dataset.xuw, dataset.yuw, dataset.numb_frames, dataset.N, dataset.t_stop)
            data["MSD"] = MSD_data["MSD"]
            data["MSD_time"] = MSD_data["MSD_time"]

        end


        if polymer_MSDbool

            polymer_MSD_data = polymer_MSD(dataset.xuw, dataset.yuw, dataset.pol_id, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop)
            data["polymer_MSD"] = polymer_MSD_data["polymer_MSD"]
            data["polymer_MSD_time"] = polymer_MSD_data["polymer_MSD_time"]

        end


        if basic_MSDbool

            basic_MSD_data = basic_MSD(dataset.xuw, dataset.yuw, dataset.numb_frames, dataset.N, dataset.t_stop)
            data["basic_MSD"] = basic_MSD_data["basic_MSD"]
            data["basic_MSD_time"] = basic_MSD_data["basic_MSD_time"]

        end


        if average_velocitybool

            average_velocity_data = average_velocity(dataset.vx, dataset.vy, dataset.numb_frames, dataset.N, dataset.t_stop)
            data["average_velocity"] = average_velocity_data["average_velocity"]
            data["average_velocity_time"] = average_velocity_data["average_velocity_time"]

        end


        if radius_of_gyrationbool
            radius_of_gyration_data = radius_of_gyration(dataset.xuw, dataset.yuw, dataset.pol_id, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop)
            data["R_2"] = radius_of_gyration_data["R_2"]
            data["R_2_time"] = radius_of_gyration_data["R_2_time"]

        end


        if end_to_end_distancebool

            end_to_end_distance_data = end_to_end_distance(dataset.xuw, dataset.yuw, dataset.pol_id, dataset.id_in_pol, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop, dataset.sizes)
            data["e_to_e_dist"] = end_to_end_distance_data["e_to_e_dist"]
            data["e_to_e_dist_time"] = end_to_end_distance_data["e_to_e_dist_time"]

        end

        close(data)

    else

        data = Dict{String, Vector{Float64}}()

        if MSDbool

            MSD_data = MSD(dataset.xuw, dataset.yuw, dataset.numb_frames, dataset.N, dataset.t_stop)
            merge!(data, MSD_data)

        end


        if polymer_MSDbool

            polymer_MSD_data = polymer_MSD(dataset.xuw, dataset.yuw, dataset.pol_id, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop)
            data["polymer_MSD"] = MSD_data["polymer_MSD"]
            data["polymer_MSD_time"] = MSD_data["polymer_MSD_time"]

        end


        if basic_MSDbool

            basic_MSD_data = basic_MSD(dataset.xuw, dataset.yuw, dataset.numb_frames, dataset.N, dataset.t_stop)
            merge!(data, basic_MSD_data)

        end


        if average_velocitybool

            average_velocity_data = average_velocity(dataset.vx, dataset.vy, dataset.numb_frames, dataset.N, dataset.t_stop)
            merge!(data, average_velocity_data)

        end


        if radius_of_gyrationbool

            radius_of_gyration_data = radius_of_gyration(dataset.xuw, dataset.yuw, dataset.pol_id, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop)
            merge!(data, radius_of_gyration_data)

        end


        if end_to_end_distancebool

            end_to_end_distance_data = end_to_end_distance(dataset.xuw, dataset.yuw, dataset.pol_id, dataset.id_in_pol, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop, dataset.sizes)
            merge!(data, end_to_end_distance_data)

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
