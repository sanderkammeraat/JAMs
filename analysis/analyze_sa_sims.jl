include("AnalysisPipeline.jl")

base_folder = base_folder = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin"





rp, sp, ap = auto_analysis_dir(base_folder, "sa_raw_data.h5"; support_raw_data_file_name_pattern = "ra_raw_data.h5")

custom_analysis_function = run_sa_analysis!

run_multithreaded_analysis(rp, ap,custom_analysis_function, support_raw_data_file_paths=sp)