include("analysisfunctions.jl")
include("loaddata.jl")



sim_folder_name = "sim_data"


ksd = [1.]
kbend = [.3]
kstretch = [1.]
fstretch = [.7]
p = [0.04, 0.06, 0.08, 0.1, 0.13, 0.15, 0.2, 0.4]
kactive = [(-1., 0.), (1., 0.), (0., 1.), (0., -1.), (1/sqrt(2), 1/sqrt(2)),(-1/sqrt(2), 1/sqrt(2)),(1/sqrt(2), -1/sqrt(2)),(-1/sqrt(2), -1/sqrt(2))]
Npol = [150]
N = [1500]

MSDbool = false
basic_MSDbool = true
average_velocitybool = true
radius_of_gyrationbool = false
end_to_end_distancebool = true


for ksd_value in ksd
for kbend_value in kbend
for kstretch_value in kstretch
for fstretch_value in fstretch
for p_value in p
for (kpar_value, kperp_value) in kactive
for Npol_value in Npol
for N_value in N 

    parameters = "p_$p_value,kpar_$kpar_value,kperp_$kperp_value" #"ksd_$ksd_value_kbend_$kbend_value_kstretch_$kstretch_value_fstretch_$f_stretch_value_p_$p_value_kperp_$kperp_value_kpar_$kpar_value_Npol_$Npol_value_N_$N_value"

    #for windows
    #path_data = joinpath("E:", "martin", sim_folder_name, parameters)

    #for linux
    path_data = joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name, parameters)

    dataset = get_data(joinpath(path_data, "raw_data.h5"))

    data = Dict{String, Vector{Float64}}()

    if MSDbool
        MSD_data = MSD(dataset.x, dataset.y, dataset.numb_frames, dataset.N, dataset.t_stop)
        merge!(data, MSD_data)
    
    if basic_MSDbool
        basic_MSD_data = basic_MSD(dataset.x, dataset.y, dataset.numb_frames, dataset.N, dataset.t_stop)
        merge!(data, basic_MSD_data)
    
    if average_velocitybool
        average_velocity_data = average_velocity(dataset.vx, dataset.vy, dataset.numb_frames, dataset.N, dataset.t_stop)
        merge!(data, average_velocity_data)
    
    if radius_of_gyrationbool
        radius_of_gyration_data = radius_of_gyration(dataset.x, dataset.y, dataset.pol_id, dataset.id_in_pol, dataset.numb_frames, dataset.N, dataset.t_stop)
        merge!(data, radius_of_gyration_data)
    
    if end_to_end_distancebool
        end_to_end_distance_data = end_to_end_distance(dataset.x, dataset.y, dataset.id_in_pol, dataset.numb_frames, dataset.Npol, dataset.N, dataset.t_stop)
        merge!(data, end_to_end_distance_data)
    end

    save_data(data, path_data)
    
end
end
end
end
end
end
end
end