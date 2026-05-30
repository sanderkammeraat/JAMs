

#base_folder = "/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin"


include("AnalysisPipeline.jl")

#base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/Nlin_20"

#base_folder = "/data2/kammeraat/sa/statistics/hex_disordered/phi_1.3/Nlin_20"

function analysis()
    #base_folder = "/data2/kammeraat/sa/statistics/hex_disordered/phi_1.3/Nlin_20"
    base_folder="/Volumes/T7_Shield/sa/statistics/hex_disordered/phi_1.3/phi_1.3/Nlin_20"
    #base_folder = "/run/media/kammeraat/T7_Shield/sa/statistics/hex_disordered/phi_1.3/vary_Nlin_partial2/"
    #base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/Nlin_20"
    rp, sp, ap = auto_analysis_dir(base_folder, "sa_raw_data.h5"; support_raw_data_file_name_pattern = "ra_raw_data.h5")


    #todo = (.!occursin.("J_0.0/", rp)) .* (.! occursin.("J_0.001", rp)) .*  (.! occursin.("J_0.01", rp)) 

    #todo = (.!occursin.("J_0.0", rp)).* (.! occursin.("J_0.12", rp)) .* (.! occursin.("J_0.16", rp)) .* (.! occursin.("J_0.1", rp)).* (.! occursin.("J_0.2", rp))

    # todo = (.!occursin.("Dr_0.01", rp)) .* (.!occursin.("J_0.0/", rp)) .* (.! occursin.("J_0.001", rp)) .*  (.! occursin.("J_0.01", rp)) 

    # rpf = rp#[todo]
    # spf = sp#[todo]
    # apf = ap#[todo]

    # for path in rpf

    #     println(path)
    # end


    # custom_analysis_function = run_sa_analysis_add_auto_p!

    #custom_analysis_function =run_sa_analysis_add_spatial_cor!

    custom_analysis_function = run_sa_analysis_add_auto_v!

    run_sequential_analysis(rp, ap,custom_analysis_function, support_raw_data_file_paths=sp,append=true)
end

analysis()

