struct Bins

    min

    max

    bin_size

    centers

    edges
    
end

struct BinnedData

    bin_centers

    N_in_bin

    bin_values

end

function create_bins(min, max, bin_size)

    edges = collect(range(min*1., step=bin_size*1., stop = max*1.))

    centers = (edges[2:end] + edges[1:end-1])/2
    return Bins(min*1., max*1., bin_size*1., centers, edges)

end

function bin_vector_data(bins::Bins, ref_vector_data, vector_data)

    N_in_bin = zeros(length(bins.centers))

    binned_vector_data = zeros(length(bins.centers))

    for i in eachindex(bins.centers)

        bool_in_bin = (ref_vector_data.> bins.edges[i]) .* (ref_vector_data.<= bins.edges[i+1])

        N_in_bin[i] = sum(bool_in_bin)

        if N_in_bin[i]>0
            binned_vector_data[i] = mean( vector_data[bool_in_bin]) 
        end
    end

    populated_bin_bool = N_in_bin .> 0

    return BinnedData(bins.centers[populated_bin_bool], N_in_bin[populated_bin_bool],binned_vector_data[populated_bin_bool])
end

#Assuming that the first dimensions corresponds to the reference data
function bin_matrix_data(bins::Bins, ref_vector_data, matrix_data)

    N_in_bin = zeros(length(bins.centers))

    binned_matrix_data = zeros((length(bins.centers),size(matrix_data)[2]))

    for i in eachindex(bins.centers)

        bool_in_bin = (ref_vector_data.> bins.edges[i]) .* (ref_vector_data.<= bins.edges[i+1])

        N_in_bin[i] = sum(bool_in_bin)

        if N_in_bin[i]>0
            binned_matrix_data[i,:] = mean( matrix_data[bool_in_bin,:], dims=1) 
        end
    end

    populated_bin_bool = N_in_bin .> 0

    return BinnedData(bins.centers[populated_bin_bool], N_in_bin[populated_bin_bool],binned_matrix_data[populated_bin_bool,:])
end






#Custom function as if having loaded a single simulation
function custom_ensemble_function_template!(ensemble_file,loaded_seed_files, seed_names)

    #Analyze stuff

    #Then return modifed analysis file
    return ensemble_file

end


#helper function

function sa_ensemble!(ensemble_file, loaded_seed_files, seed_names)

    #Use the first seed to extract system information

    reference_seed = loaded_seed_files[1]

    ensemble_file["Nseeds"] = length(loaded_seed_files)

    ensemble_file["seed_names"] = seed_names

    ensemble_file["t"] = reference_seed["t"]

    ensemble_file["dt"] = reference_seed["dt"]

    ensemble_file["Nt"] = reference_seed["Nt"]

    ensemble_file["v0"] = reference_seed["v0"]

    ensemble_file["Dr"] = reference_seed["Dr"]

    ensemble_file["J"] = reference_seed["J"]

    ensemble_file["k"] = reference_seed["k"]

    ensemble_file["Nint"] = reference_seed["Nint"]

    min_t_ind = reference_seed["min_t_ind"]
    ensemble_file["min_t_ind"] = reference_seed["min_t_ind"]

    ensemble_file["type"] = reference_seed["type"]

    ensemble_file["R"] = reference_seed["R"]



    create_group(ensemble_file, "eigenmodes")
    create_group(ensemble_file["eigenmodes"],"eigvals" )
    for i in eachindex(loaded_seed_files)
        ensemble_file["eigenmodes"]["eigvals"][seed_names[i]]=loaded_seed_files[i]["eigenmodes"]["eigvals"]
    end

    eigvalbins = create_bins(0, 1,0.002)


    ensemble_file["eigval_bin_centers"] = eigvalbins.centers

    v_projs = []
    eigvals = []
    for i in eachindex(loaded_seed_files)
        v_projs=vcat(v_projs, mean(reference_seed["projs"]["v_projs"][:,min_t_ind:end].^2, dims=2)[:,1])
        eigvals=vcat(eigvals, loaded_seed_files[i]["eigenmodes"]["eigvals"])
    end

    binned_v_projs_data = bin_vector_data(eigvalbins, eigvals, v_projs)
    create_group(ensemble_file, "v_projs_time_avg")

    ensemble_file["v_projs_time_avg"]["eigval_bin_centers"] = binned_v_projs_data.bin_centers
    ensemble_file["v_projs_time_avg"]["v_projs_time_avg"] = binned_v_projs_data.bin_values
    ensemble_file["v_projs_time_avg"]["N_in_bin"] = binned_v_projs_data.N_in_bin


    vrms = []
    for i in eachindex(loaded_seed_files)
        vrms=vcat(vrms, loaded_seed_files[i]["vrms_particle_avg_time_avg"])
    end
    ensemble_file["vrms"]= mean(vrms)


    rbins = create_bins(0, 50,2.5)
    ensemble_file["r_bin_centers"] = rbins.centers

    vrms_r = []
    r = []
    for i in eachindex(loaded_seed_files)
        r=vcat(r, loaded_seed_files[i]["r_particle_time_avg"])
        vrms_r = vcat(vrms_r,loaded_seed_files[i]["vrms_particle_time_avg"] )
    end
    create_group(ensemble_file, "vrms_r")

    binned_vrms_r_data = bin_vector_data(rbins, r, vrms_r)
    ensemble_file["vrms_r"]["r_bin_centers"]= binned_vrms_r_data.bin_centers
    ensemble_file["vrms_r"]["vrms_r"]= binned_vrms_r_data.bin_values
    ensemble_file["vrms_r"]["N_in_bin"]= binned_vrms_r_data.N_in_bin

    

    X = []
    w = []
    for i in eachindex(loaded_seed_files)

        if i==1

            w =loaded_seed_files[i]["FT_v_projs"]["w"]

            X = loaded_seed_files[i]["FT_v_projs"]["Xf2"]
        else

            X = vcat(X, loaded_seed_files[i]["FT_v_projs"]["Xf2"])
        end
    end

    binned_X = bin_matrix_data(eigvalbins, eigvals, X)

    #display("Collecting FTs")
    create_group(ensemble_file, "FT_v_projs")

    ensemble_file["FT_v_projs"]["w"] = w
    ensemble_file["FT_v_projs"]["eigval_bin_centers"] = binned_X.bin_centers
    ensemble_file["FT_v_projs"]["X2"] = binned_X.bin_values

    #auto_p
    Cavg = []
    deltat = []
    for i in eachindex(loaded_seed_files)

        if i==1

            deltat =loaded_seed_files[i]["auto_p"]["deltat"]

            Cavg = loaded_seed_files[i]["auto_p"]["Cavg"]
        else

            Cavg = vcat(Cavg, loaded_seed_files[i]["auto_p"]["Cavg"])
        end
    end
    create_group(ensemble_file, "auto_p")
    ensemble_file["auto_p"]["Cavg"] = mean(Cavg, dims=1)[1,:] 
    ensemble_file["auto_p"]["deltat"] = deltat


    return ensemble_file
end

function sa_ensemble_add_FT_px!(ensemble_file, loaded_seed_files, seed_names)


    X = []
    w = []
    for i in eachindex(loaded_seed_files)

        if i==1

            w =loaded_seed_files[i]["FT_px"]["w"]

            X = reshape(loaded_seed_files[i]["FT_px"]["pavg_X2"], 1, length(w))
        else

            X = vcat(X, reshape(loaded_seed_files[i]["FT_px"]["pavg_X2"], 1, length(w)))
        end
    end
    create_group(ensemble_file,"FT_px_w")
    ensemble_file["FT_px_w"]["X2"] = mean(X, dims=1)[1,:]
    ensemble_file["FT_px_w"]["w"] = w


    return ensemble_file

end

