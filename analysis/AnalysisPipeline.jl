using Distributed

include(joinpath("..","src","Engine.jl"))
include("AnalysisFunctions.jl")
include("CustomAnalysisFunctions.jl")
include("CustomEnsembleFunctions.jl")

#Helper function
function load_file(file_location)

    file = jldopen(file_location, "r",iotype=IOStream)

end

#Helper function
function initialize_file(save_path; overwrite = false, append=false)

    save_folder_path = mkpath(dirname(save_path))

    if append && overwrite
        error("Set only append or error to true, not both!")
    else
        if !isfile(save_path)

            file = h5open(save_path,"cw")

            return file

        elseif isfile(save_path) && append
            file = h5open(save_path,"cw")

            return file

        elseif isfile(save_path) && overwrite
            file = h5open(save_path,"w")

            return file

        else
            error("Initialized file already existing, aborting to prevent overwriting.")
        end
    end
end


#With data_dict, we mean a Julia dict in which we will collect data that we need for analysis
#With analysis_file we mean an initialized analysis file in which we will store the analysis results
#The idea is to first expose the necessary data from the simulation file to the datadict, and then let the add_analysis files use this to calculate something.

#General analysis
function analyze_single(raw_data_file_path, analysis_save_path, custom_analysis_function; support_raw_data_file_path=nothing,  overwrite = false, append = false)


    #Preparing files
    loaded_raw_data_file = load_file(raw_data_file_path)

    loaded_support_raw_data_file = !isnothing(support_raw_data_file_path) ?  load_file(support_raw_data_file_path) : nothing

    try 
        initialized_analysis_file = initialize_file(analysis_save_path; overwrite = overwrite, append=append)

        try 
            
            #Actually running the analysis
            initialized_analysis_file = custom_analysis_function(initialized_analysis_file,loaded_raw_data_file; support_raw_data_file = loaded_support_raw_data_file)

            #Closing the files
            close(initialized_analysis_file)

            if !isnothing(loaded_support_raw_data_file)
                close(loaded_support_raw_data_file)
            end

            close(loaded_raw_data_file)

        catch e 

            println("Something in the custom analysis function went wrong. Safely aborting and closing the files.")
            close(initialized_analysis_file)

            if !isnothing(loaded_support_raw_data_file)
                close(loaded_support_raw_data_file)
            end

            close(loaded_raw_data_file)

            rethrow(e)
        end
    catch e
        println("Something went wrong. Safely aborting and closing the files.")

        if !isnothing(loaded_support_raw_data_file)
            close(loaded_support_raw_data_file)
        end

        close(loaded_raw_data_file)

        rethrow(e)

    end

end

function ensemble_single(seed_file_paths, custom_ensemble_function; overwrite = false, append = false)


    #Preparing files

    loaded_seed_files = [ load_file(seed_file_path) for seed_file_path in seed_file_paths]

    seed_names = [ splitpath(seed_file_path)[end] for seed_file_path in seed_file_paths ]

    try 
        initialized_ensemble_file = initialize_file(ensemble_save_path; overwrite = overwrite, append=append)

        try 
            
            #Actually running the analysis
            initialized_ensemble_file = custom_ensemble_function(initialized_ensemble_file,loaded_seed_files, seed_names)

            #Closing the files
            close(initialized_ensemble_file)

            close.(loaded_seed_files)

        catch e 

            println("Something in the custom analysis function went wrong. Safely aborting and closing the files.")
            close(initialized_ensemble_file)

            close.(loaded_seed_files)

            rethrow(e)
        end
    catch e
        println("Something went wrong. Safely aborting and closing the files.")


        close.(loaded_seed_files)

        rethrow(e)

    end

end


function movie_single(raw_data_file_path, movie_save_path, custom_movie_function)


    #Preparing files
    loaded_raw_data_file = load_file(raw_data_file_path)

    try 
        mkpath(dirname(movie_save_path))
        custom_movie_function(loaded_raw_data_file,movie_save_path)


        close(loaded_raw_data_file)

    catch e
        println("Something went wrong. Safely aborting and closing the raw data file.")

        close(loaded_raw_data_file)

        rethrow(e)

    end

end



#helper function
function save_dict2h5!(file,dict, dictsavename)

    create_group(file, dictsavename)

    for (key, val) in dict
        file[dictsavename][string(key)] = val
    end
    return file
end

#helper function
function extract_frame_data_for_type(datakey, type, frame_data)

    return frame_data[datakey][ frame_data["type"].==type ]

end

function findfile(directory, filepattern)

    paths = String[]

    for (root, dirs, files) in walkdir(directory)

        for file in files

            if occursin(filepattern, file) && occursin(".h5", file)

                push!(paths, joinpath(root, file))
            end
        end
    end

    return paths
end



function auto_analysis_dir(base_folder, raw_data_file_name_pattern; support_raw_data_file_name_pattern = nothing, mimic_source_dir_name = "simdata",mimic_target_dir_name = "analysis")

    raw_data_file_paths = findfile(base_folder, raw_data_file_name_pattern)

    support_raw_data_file_paths = !isnothing(support_raw_data_file_name_pattern) ? findfile(base_folder, support_raw_data_file_name_pattern) : nothing

    seed_dir_names = [ splitpath(dirname(file_path))[end] for file_path in raw_data_file_paths]

    analysis_save_folders = replace.( dirname.(dirname.(raw_data_file_paths)), mimic_source_dir_name => mimic_target_dir_name) 

    analysis_save_paths = [joinpath(analysis_save_folders[i],seed_dir_names[i]*".h5") for i=eachindex(seed_dir_names)]


    return raw_data_file_paths, support_raw_data_file_paths, analysis_save_paths

end

function auto_ensemble_dir(base_folder, analysis_file_name_pattern; mimic_source_dir_name = "analysis", mimic_target_dir_name = "ensembles")

    #Find all analysis files
    analysis_file_paths = findfile(base_folder, analysis_file_name_pattern)

    display(analysis_file_paths)

    #Construct ensemble save folder structures
    ensemble_save_folders = unique(replace.( dirname.(analysis_file_paths), mimic_source_dir_name => mimic_target_dir_name))

    ensemble_save_paths = joinpath.(ensemble_save_folders,"ensemble.h5")

    seed_file_paths_per_ensemble_folder = []
    
    for i in eachindex(ensemble_save_folders)

        #Corresponding analysis folder
        analysis_folder  = replace.( ensemble_save_folders[i], mimic_target_dir_name => mimic_source_dir_name) 

        #Find all seeds in this folder
        seed_paths_in_ensemble = findfile(analysis_folder, analysis_file_name_pattern)

        #Add these file paths to the corresponding ensemble
        push!(seed_file_paths_per_ensemble_folder,seed_paths_in_ensemble)
    end


    return ensemble_save_paths, seed_file_paths_per_ensemble_folder
end

function auto_movie_dir(base_folder, raw_data_file_name_pattern, mimic_source_dir_name = "simdata", mimic_target_dir_name = "movies")

    raw_data_file_paths = findfile(base_folder, raw_data_file_name_pattern)

    seed_dir_names = [ splitpath(dirname(file_path))[end] for file_path in raw_data_file_paths]

    movie_save_folders = replace.( dirname.(dirname.(raw_data_file_paths)), mimic_source_dir_name => mimic_target_dir_name) 

    movie_save_paths = [joinpath(movie_save_folders[i],seed_dir_names[i]*".mp4") for i=eachindex(seed_dir_names)]


    return raw_data_file_paths, movie_save_paths

end


function run_distributed_analysis(raw_data_file_paths, analysis_save_paths,custom_analysis_function; support_raw_data_file_paths=nothing,  overwrite = false, append = false)

    @sync @distributed for i in eachindex(raw_data_file_paths)
        
        raw_data_file_path = raw_data_file_paths[i]
        println(raw_data_file_path)
        analysis_save_path = analysis_save_paths[i]
        support_raw_data_file_path = !isnothing(support_raw_data_file_paths) ? support_raw_data_file_paths[i] : nothing

        analyze_single(raw_data_file_path, analysis_save_path, custom_analysis_function; support_raw_data_file_path=support_raw_data_file_path,  overwrite = overwrite, append = append)

    end

end

function run_multithreaded_analysis(raw_data_file_paths, analysis_save_paths,custom_analysis_function; support_raw_data_file_paths=nothing,  overwrite = false, append = false)

    Threads.@threads for i in eachindex(raw_data_file_paths)
        raw_data_file_path = raw_data_file_paths[i]
        println(raw_data_file_path)
        analysis_save_path = analysis_save_paths[i]
        support_raw_data_file_path = !isnothing(support_raw_data_file_paths) ? support_raw_data_file_paths[i] : nothing

        analyze_single(raw_data_file_path, analysis_save_path, custom_analysis_function; support_raw_data_file_path=support_raw_data_file_path,  overwrite = overwrite, append = append)

    end

end

function run_sequential_analysis(raw_data_file_paths, analysis_save_paths,custom_analysis_function; support_raw_data_file_paths=nothing,  overwrite = false, append = false)

    for i in eachindex(raw_data_file_paths)

        raw_data_file_path = raw_data_file_paths[i]
        println(raw_data_file_path)
        analysis_save_path = analysis_save_paths[i]
        support_raw_data_file_path = !isnothing(support_raw_data_file_paths) ? support_raw_data_file_paths[i] : nothing

        analyze_single(raw_data_file_path, analysis_save_path, custom_analysis_function; support_raw_data_file_path=support_raw_data_file_path,  overwrite = overwrite, append = append)

    end

end

function run_multithreaded_movie(raw_data_file_paths, movie_save_paths,custom_movie_function)

    Threads.@threads for i in eachindex(raw_data_file_paths)
        raw_data_file_path = raw_data_file_paths[i]
        println(raw_data_file_path)
        movie_save_path = movie_save_paths[i]

        movie_single(raw_data_file_path, movie_save_path, custom_movie_function)
    end

end

function run_sequential_movie(raw_data_file_paths, movie_save_paths,custom_movie_function; only_seed="seed_1")

    @showprogress showspeed=true for i in eachindex(raw_data_file_paths)

        raw_data_file_path = raw_data_file_paths[i]

        if splitpath(raw_data_file_path)[end-1] == only_seed
        
            println(raw_data_file_path)
            movie_save_path = movie_save_paths[i]

            movie_single(raw_data_file_path, movie_save_path, custom_movie_function)
        end
    end

end
