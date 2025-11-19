#Custom function as if having loaded a single simulation
function custom_ensemble_function(ensemble_file,loaded_seed_files, seed_names)

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

    

    ensemble_file["min_t_ind"] = reference_seed["min_t_ind"]

    ensemble_file["type"] = reference_seed["type"]

    ensemble_file["R"] = reference_seed["R"]


    create_group(ensemble_file, "v_projs")
    create_group(ensemble_file["v_projs"], "seeds")

    for seed_name in seed_names

        create_group(ensemble_file["v_projs"]["seeds"], seed_name)
        #save_dict2h5!(ensemble_file["v_projs"]["seeds"][seed_name],v_proj,)

    end




    return ensemble_file
end