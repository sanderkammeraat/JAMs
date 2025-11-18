


include("AnalysisPipeline.jl")



raw_data_file_path = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin/Nlin_10/simdata/Dr_0.01/J_0.1/seed_1/sa_raw_data.h5"

support_raw_data_file_path = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin/Nlin_10/simdata/Dr_0.01/J_0.1/seed_1/ra_raw_data.h5"

analysis_save_path = joinpath(homedir(), "test_new_analysis", "test.h5")


analyze_single(raw_data_file_path, analysis_save_path, run_sa_analysis!; support_raw_data_file_path=support_raw_data_file_path, overwrite=true)


base_folder = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin"



raw_data_file_paths, support_raw_data_file_paths, analysis_save_paths = auto_analysis_dir(base_folder, "sa_raw_data.h5", support_raw_data_file_name_pattern = "ra_raw_data.h5" )



readdir(base_folder)

walkdir(base_folder)


# analysis = load_file(analysis_save_path)

# analysis["projs"]["v_projs"]

# vrms=  analysis["vrms_particle_time_avg"]


# r = analysis["r_particle_time_avg"]

# plot(r, vrms/analysis["v0"])
# ylims!(0,1)

# close(analysis)