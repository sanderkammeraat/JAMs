include("analysisfunctions.jl")
include("loaddata.jl")



sim_folder_name = "sim_data"


ksd = [1.]
kbend = [.3]
kstretch = [1.]
fstretch = [.7]
p = [0.01, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3]
kperp = [0]
kpar = [-1]
Npol = [150]
N = [1500]

MSDbool = true
average_velocitybool = true
radius_of_gyrationbool = false


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
    #path_raw_data = joinpath("E:", "martin", sim_folder_name, parameters, "raw_data.h5")

    #for linux
    path_raw_data = joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name, parameters, "raw_data.h5")

    data = get_data(path_raw_data)

    if MSDbool
        save_data(MSD(data.x, data.y, data.numb_frames, data.N), "MSD", joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name, parameters))
    end
    if average_velocitybool
        save_data(average_velocity(data.vx, data.vy, data.numb_frames, data.N, convert(Int64, floor(data.numb_frames/100))+1), "average_velocity", joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name, parameters))
    end
    if radius_of_gyrationbool
        save_data(radius_of_gyration(data.x, data.y, data.pol_id, data.id_in_pol, data.numb_frames, data.N), "radius_of_gyration", joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name, parameters))

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

