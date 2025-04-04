

function readdir_filt(folder_path)
    directories = filter!(e->e!=".DS_Store",readdir(folder_path))
    return directories

end

function construct_folder_tree_param_param_seed(base_folder)
    tree = Dict()

    #Base folders (first parameter axis)
    dirs = readdir_filt(base_folder)
    
    #Subfolder (second parameter axis)
    for dir in dirs
        tree[dir]=Dict()
        subdirs = readdir_filt(joinpath( base_folder, dir ))
    
        for subdir in  subdirs
            tree[dir][subdir]=Dict()
    
            seeddirs =  readdir_filt(joinpath( base_folder, dir, subdir ))
    
            for seeddir in seeddirs
                tree[dir][subdir][seeddir]= joinpath( base_folder, dir, subdir, seeddir )
            end
        end
        
    
    end
    return tree
end

function initialize_analysis_file!(raw_data_file,analysis_file)

    integration_info = raw_data_file["integration_info"]

    system = read(raw_data_file,"system")

    for info in keys(integration_info)
        analysis_file["integration_info/$(info)"] = raw_data_file["integration_info"][info]
    end

    for forcetype in keys(system["forces"])

        for force in keys(system["forces"][forcetype])

            for force_param in keys(system["forces"][forcetype][force])

                analysis_file["system/forces/$(forcetype)/$(force)/$(force_param)"]=system["forces"][forcetype][force][force_param]
            end
        end
    end

    for dofevolver in keys(system["dofevolvers"])
        analysis_file["system/dofevolvers/$(dofevolver)"] = system["dofevolvers"][dofevolver]
    end

end
#Does nothing else than opening and saving analysis file.
# analyze_single_seed_inner is a function that should take  (analysis_file, system, integration_info, frames) as arguments and return analysis_file.
function analyze_single_seed_outer(raw_data_file_path, analysis_file_path, analyze_single_seed_inner; overwrite=false)

    raw_data_file = jldopen(raw_data_file_path, "r")

    if overwrite
        analysis_file = jldopen(analysis_file_path, "w")
        initialize_analysis_file!(raw_data_file,analysis_file)

    else
        if !isfile(analysis_file_path)
            
            jldopen(analysis_file_path, "a+") do analysis_file
                initialize_analysis_file!(raw_data_file,analysis_file)
            end
        end
        #But for the analysis data, so reload in append only mode 
        analysis_file = jldopen(analysis_file_path, "a+")
    end

    #Expose by default the following data
    system = analysis_file["system"]
    integration_info  = analysis_file["integration_info"]
    frames = raw_data_file["frames"]

    #Do stuff with custom inner analysis function

    try
        analysis_file = analyze_single_seed_inner(analysis_file, system, integration_info, frames)
    catch e
        close(analysis_file)
        close(raw_data_file)
        error("JAMS: data(group) already present in analysis file, aborted to prevent overwriting.  Suggested solution: check if overwriting is desired and if so, use overwrite=true.)")
        rethrow(e)

    end

    #Now close
    close(analysis_file)
    close(raw_data_file)
end

function run_serial_analysis_param1_param2_seed(tree, analyze_single_seed_inner, analysis_base_folder; raw_data_file_name="raw_data.jld2", overwrite=false)

     for (param1, subdict) in tree

        for (param2, seeddict) in subdict

            @showprogress dt = 1 desc="Analysis in progress..." showspeed=true for (seed, seedpath) in  seeddict

                raw_data_file_path = joinpath(seedpath, raw_data_file_name)

                analysis_file_path = joinpath( mkpath(joinpath(analysis_base_folder,param1, param2)), "$(seed).jld2")

                analyze_single_seed_outer(raw_data_file_path,analysis_file_path, analyze_single_seed_inner, overwrite=overwrite)
            end

        end

    end
end

function run_multithreaded_analysis_param1_param2_seed(tree, analyze_single_seed_inner, analysis_base_folder; raw_data_file_name="raw_data.jld2", overwrite=false)

    k1s = collect(keys(tree))
    for k1 in k1s
        k2s = collect(keys(tree[k1]))

        Threads.@threads for k2 in k2s
        k3s = collect(keys(tree[k1][k2]))

           @showprogress dt = 1 desc="Analysis in progress..." showspeed=true for k3 in  k3s
            
               seedpath = tree[k1][k2][k3]
               seed = k3
               raw_data_file_path = joinpath(seedpath, raw_data_file_name)

               analysis_file_path = joinpath( mkpath(joinpath(analysis_base_folder,k1, k2)), "$(seed).jld2")

               analyze_single_seed_outer(raw_data_file_path,analysis_file_path, analyze_single_seed_inner, overwrite=overwrite)
           end

       end

   end
end


function extract_frame_data_for_type(datakey, type, frame_data)

    return frame_data[datakey][ frame_data["type"].==type ]

end

function acces_param1_param2_seedanalysis(tree, dofunctions)

    for (param1, subdict) in tree

        for (param2, seeddict) in subdict

            for (seed, seedpath) in  seeddict

                jldopen(seedpath, "r") do seedanalysis_file

                    for dofunction in dofunctions
                        dofunction(seed, seedanalysis_file)
                    end

                end

            end
        end
    end

end
