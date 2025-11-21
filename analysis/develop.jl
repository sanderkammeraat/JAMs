



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

seed_file_paths = sf[3]

loaded_seed_files = [ load_file(seed_file_path) for seed_file_path in seed_file_paths]


reference_seed = loaded_seed_files[3]

min_t_ind = reference_seed["min_t_ind"]

v_projs_time_avg =  mean(reference_seed["projs"]["v_projs"][:,min_t_ind:end].^2, dims=2)[:,1]

eigvals = reference_seed["eigenmodes"]["eigvals"]


plot(eigvals, log10.(v_projs_time_avg))


#close.(loaded_seed_files)








bins = create_bins(0, 10,0.002)
bins.centers
bins.edges



f = Figure()

ax = Axis(f[1,1], xlabel=L"w", ylabel=L"FT", yscale=log10, xscale=log10)
min_t_ind = 500
X = []
eigvals = []
w = []
for i in eachindex(loaded_seed_files)

    if i==1
        w =vcat(w, loaded_seed_files[i]["FT_v_projs"]["w"])

        
        X = loaded_seed_files[i]["FT_v_projs"]["Xf2"]


        display(size(X))
        display(size(loaded_seed_files[i]["eigenmodes"]["eigvals"]))
        eigvals=vcat(eigvals, loaded_seed_files[i]["eigenmodes"]["eigvals"])
    else
    eigvals=vcat(eigvals, loaded_seed_files[i]["eigenmodes"]["eigvals"])

        scatter!(ax,loaded_seed_files[i]["FT_v_projs"]["Xf2"][40,:])
        display(f)
    

    X = vcat(X, loaded_seed_files[i]["FT_v_projs"]["Xf2"])
    end
end

eigvals
X

binned_X = bin_matrix_data(bins, eigvals, X)


using GLMakie
f = Figure()
ax = Axis(f[1,1], xlabel=L"w", ylabel=L"FT", yscale=log10)
for i = 1:10

    lines!(ax, w, binned_X.bin_values[i,:], color = i, colorrange=(0,10))
end
display(f)

image(log10.(transpose(binned_X.bin_values)), colormap=:viridis)



A = [1 2 ; 3 4]

A = vcat(A, [5 6 ; 7 8])

A[A[:,1].>1,:]


binned_data = bin_vector_data(bins, eigvals, v_projs)


plot(binned_data.bin_centers, binned_data.bin_values)

close.(loaded_seed_files)


#%%

