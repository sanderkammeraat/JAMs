

#base_folder = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin"


include("AnalysisPipeline.jl")

#base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/Nlin_20"
#
#base_folder = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/fa_extra/"
#
#base_folder = "/Volumes/T7_Shield/sa/statistics/free/phi_1.3/N_2000/"

base_folder = "/Volumes/T7_Shield/sa/statistics/free/phi_1.3/fa"


#base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/vary_Nlin"


ef, sf = auto_ensemble_dir(base_folder,"seed_")


#base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/Nlin_20"


run_sequential_ensemble(ef, sf, sa_ensemble_free!, overwrite=false)



