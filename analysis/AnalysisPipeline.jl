using Distributed
include(joinpath("..","src","Engine.jl"))
include("AnalysisFunctions.jl")
include("CustomAnalysisFunctions.jl")

#Helper function
function load_file(file_location)

    file = jldopen(file_location, "r",iotype=IOStream)

end

#Helper function
function initialize_analysis_file(save_path; overwrite = false, append=false)

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
            error("Analysis file already existing, aborting to prevent overwriting.")
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
        initialized_analysis_file = initialize_analysis_file(analysis_save_path; overwrite = overwrite, append=append)

        try 
            
            #Actually running the analysis
            initialized_analysis_file = custom_analysis_function(initialized_analysis_file,loaded_raw_data_file; support_raw_data_file = loaded_support_raw_data_file)

            #Closing the files
            close(initialized_analysis_file)

            close(loaded_support_raw_data_file)

            close(loaded_raw_data_file)

        catch e 

            println("Something in the custom analysis function went wrong. Safely aborting and closing the files.")
            close(initialized_analysis_file)

            close(loaded_support_raw_data_file)

            close(loaded_raw_data_file)

            rethrow(e)
        end
    catch e
        println("Something went wrong. Safely aborting and closing the files.")

        close(loaded_support_raw_data_file)

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

#helper function, thanks to https://discourse.julialang.org/t/find-all-files-named-findthis-csv-in-nested-subfolders-of-rootfolder/118096/7
findfile(directory, file) = [joinpath(root, file) for (root, dirs, files) in walkdir(directory) if file in files]
#

function auto_analysis_dir(base_folder, raw_data_file_name_pattern; support_raw_data_file_name_pattern = nothing, mimic_source_dir_name = "simdata",mimic_target_dir_name = "analysis")

    raw_data_file_paths = findfile(base_folder, raw_data_file_name_pattern)

    support_raw_data_file_paths = !isnothing(support_raw_data_file_name_pattern) ? findfile(base_folder, support_raw_data_file_name_pattern) : nothing

    seed_dir_names = [ splitpath(dirname(file_path))[end] for file_path in raw_data_file_paths]

    analysis_save_folders = replace.( dirname.(dirname.(raw_data_file_paths)), mimic_source_dir_name => mimic_target_dir_name) 

    analysis_save_paths = [joinpath(analysis_save_folders[i],seed_dir_names[i]*".h5") for i=eachindex(seed_dir_names)]


    return raw_data_file_paths, support_raw_data_file_paths, analysis_save_paths

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

