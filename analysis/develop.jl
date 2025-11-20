



include("AnalysisPipeline.jl")



# raw_data_file_path = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin/Nlin_10/simdata/Dr_0.01/J_0.1/seed_1/sa_raw_data.h5"

# support_raw_data_file_path = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin/Nlin_10/simdata/Dr_0.01/J_0.1/seed_1/ra_raw_data.h5"

# analysis_save_path = joinpath(homedir(), "test_new_analysis", "test.h5")


# analyze_single(raw_data_file_path, analysis_save_path, run_sa_analysis!; support_raw_data_file_path=support_raw_data_file_path, overwrite=true)


# base_folder = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin"




# raw_data_file_paths, support_raw_data_file_paths, analysis_save_paths = auto_analysis_dir(base_folder, "sa_raw_data.h5", support_raw_data_file_name_pattern = "ra_raw_data.h5" )


base_folder = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin"


ef, sf = auto_ensemble_dir(base_folder,"seed_")


ef[2]

sf[2][1]

seed_file_paths = sf[2]

loaded_seed_files = [ load_file(seed_file_path) for seed_file_path in seed_file_paths]


reference_seed = loaded_seed_files[1]

min_t_ind = reference_seed["min_t_ind"]

v_projs_time_avg =  mean(reference_seed["projs"]["v_projs"][:,min_t_ind:end].^2, dims=2)[:,1]

eigvals = reference_seed["eigenmodes"]["eigvals"]


plot(log10.(eigvals), log10.(v_projs_time_avg))


#close.(loaded_seed_files)








bins = create_bins(0, 10,0.05)
bins.centers
bins.edges





v_projs = []
eigvals = []
for i in eachindex(loaded_seed_files)

    v_projs=vcat(v_projs, mean(reference_seed["projs"]["v_projs"][:,min_t_ind:end].^2, dims=2)[:,1])
    eigvals=vcat(eigvals, loaded_seed_files[i]["eigenmodes"]["eigvals"])
end




binned_data = bin_vector_data(bins, eigvals, v_projs)


plot(binned_data.bin_centers, binned_data.bin_values)

close.(loaded_seed_files)


#%%

