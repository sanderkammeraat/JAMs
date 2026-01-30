

#base_folder = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin"


include("AnalysisPipeline.jl")

#base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/Nlin_20"

base_folder = "/data2/kammeraat/sa/statistics/hex_disordered/phi_1.3/Nlin_20"

ef, sf = auto_ensemble_dir(base_folder,"seed_")


#base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/Nlin_20"


run_sequential_ensemble(ef, sf, sa_ensemble!, overwrite=false)



