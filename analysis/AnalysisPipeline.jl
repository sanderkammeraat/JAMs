include(joinpath("..","src","Engine.jl"))
include("AnalysisFunctions.jl")

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

#Custom function as if having loaded a single simulation
function run_custom_analysis_function_template!(analysis_file, raw_data_file; support_raw_data_file = nothing)

    #Analyze stuff

    #Then return modifed analysis file
    return analysis_file

end



function run_sa_analysis!(analysis_file, raw_data_file; support_raw_data_file = nothing)


    frames = raw_data_file["frames"]

    system = raw_data_file["system"]

    frames_support = support_raw_data_file["frames"]

    t =  raw_data_file["integration_info"]["save_tax"]
    analysis_file["t"] = t

    v0 = frames["1"]["v0"][1]
    analysis_file["v0"] = v0

    Dr = frames["1"]["Dr"][1]
    analysis_file["Dr"] = Dr

    J = system["forces"]["external"]["self_align_with_v_unit_force"]["β"]
    analysis_file["J"] = J


    k = system["forces"]["pair"]["soft_disk_force"]["karray"]
    analysis_file["k"] = k

    R = frames["1"]["R"]
    analysis_file["R"] = R


    type = frames["1"]["type"]
    analysis_file["type"] = type


    Nt = length(t)
    analysis_file["Nt"] = Nt
    #Ignore boundary particles
    Nint = length(extract_frame_data_for_type("id",1,frames["1"]))
    analysis_file["Nint"] = Nint

    min_t_ind = 500
    analysis_file["min_t_ind"] = min_t_ind

    dt = t[2] - t[1]
    analysis_file["dt"] = dt


    x = zeros(Nint, Nt)
    y = zeros(Nint, Nt)

    vx = zeros(Nint, Nt)
    vy = zeros(Nint, Nt)

    px = zeros(Nint, Nt)
    py = zeros(Nint, Nt)

    qx = zeros(Nint, Nt)
    qy = zeros(Nint, Nt)


    @views for i in 1:Nt
        x[:,i] .= extract_frame_data_for_type("x", 1, frames[string(i)])
        y[:,i] .= extract_frame_data_for_type("y", 1, frames[string(i)])

        vx[:,i] .= extract_frame_data_for_type("vx", 1, frames[string(i)])
        vy[:,i] .= extract_frame_data_for_type("vy", 1, frames[string(i)])

        px[:,i] .= extract_frame_data_for_type("px", 1, frames[string(i)])
        py[:,i] .= extract_frame_data_for_type("py", 1, frames[string(i)])

        qx[:,i] .= extract_frame_data_for_type("qx", 1, frames[string(i)])
        qy[:,i] .= extract_frame_data_for_type("qy", 1, frames[string(i)])

    end


    x0 = frames_support[string(length(frames_support))]["x"]
    y0 = frames_support[string(length(frames_support))]["y"]


    #Find out oscillation frequency of particles from p
    FT_px_pp = temporal_Fourier_transform(dt, px, min_t_ind = min_t_ind, output_not_avg=true)
    FT_py_pp = temporal_Fourier_transform(dt, py, min_t_ind = min_t_ind, output_not_avg=true)
    create_group(analysis_file, "FT_p_per_particle")
    analysis_file["FT_p_per_particle"]["FT_px_w_max"] = FT_px_pp["w_max"]
    analysis_file["FT_p_per_particle"]["FT_py_w_max"] = FT_py_pp["w_max"]


    #Find out oscillation frequency of particles from v
    FT_vx_pp = temporal_Fourier_transform(dt, vx, min_t_ind = min_t_ind, output_not_avg=true)
    FT_vy_pp = temporal_Fourier_transform(dt, vy, min_t_ind = min_t_ind, output_not_avg=true)
    create_group(analysis_file, "FT_v_per_particle")
    analysis_file["FT_v_per_particle"]["FT_vx_w_max"] = FT_vx_pp["w_max"]
    analysis_file["FT_v_per_particle"]["FT_vy_w_max"] = FT_vy_pp["w_max"]

    #Calculate particle average spectrum of p
    FT_px = temporal_Fourier_transform(dt, px, min_t_ind = min_t_ind)
    FT_py = temporal_Fourier_transform(dt, py, min_t_ind = min_t_ind)
    save_dict2h5!(analysis_file,FT_px, "FT_px")
    save_dict2h5!(analysis_file,FT_py, "FT_py")

    #Calculate particle average spectrum of v
    FT_vx = temporal_Fourier_transform(dt, vx, min_t_ind = min_t_ind)
    FT_vy = temporal_Fourier_transform(dt, vy, min_t_ind = min_t_ind)
    save_dict2h5!(analysis_file,FT_vx, "FT_vx")
    save_dict2h5!(analysis_file,FT_vy, "FT_vy")

    #Calculate mean radial distance of particles
    r_particle_time_avg = mean(sqrt.(x.^2 + y.^2)[:,min_t_ind:end],dims=2)[:,1]
    analysis_file["r_particle_time_avg"] = r_particle_time_avg

    #vrms info 
    vrms_particle_time_avg = mean(sqrt.(vx.^2 + vy.^2)[:,min_t_ind:end],dims=2)[:,1]
    analysis_file["vrms_particle_time_avg"] = vrms_particle_time_avg

    vrms_particle_avg_time_avg = mean( sqrt.(vx.^2 + vy.^2)[:,min_t_ind:end] )
    analysis_file["vrms_particle_avg_time_avg"] =  vrms_particle_avg_time_avg


    vrms_particle_avg_time = mean(sqrt.(vx.^2 + vy.^2)[:,min_t_ind:end],dims=1)[1,:]
    analysis_file["vrms_particle_avg_time"] = vrms_particle_avg_time


    #avg p info
    px_particle_time_avg = mean(px[:,min_t_ind:end],dims=2)[:,1]
    analysis_file["px_particle_time_avg"] = px_particle_time_avg

    py_particle_time_avg = mean(py[:,min_t_ind:end],dims=2)[:,1]
    analysis_file["py_particle_time_avg"] = py_particle_time_avg

    px_particle_avg_time_avg = mean( px[:,min_t_ind:end] )
    analysis_file["px_particle_avg_time_avg"] =  px_particle_avg_time_avg

    py_particle_avg_time_avg = mean( py[:,min_t_ind:end] )
    analysis_file["py_particle_avg_time_avg"] =  py_particle_avg_time_avg


    px_particle_avg_time = mean(px[:,min_t_ind:end],dims=1)[1,:]
    analysis_file["px_particle_avg_time"] = px_particle_avg_time

    py_particle_avg_time = mean(py[:,min_t_ind:end],dims=1)[1,:]
    analysis_file["py_particle_avg_time"] = py_particle_avg_time


    #Heavy work: doing the mode analysis
    #With boundary
    display("Constructing D_wb")
    D_wb = construct_D(x0, y0, k, R, type)


    #@profview_allocs D_wb = construct_D(x0, y0, k, R, type)
    #Boundary particles are stored at the end
    interior_indmax = sum(type.==1)

    display("Slicing out D")
    D = Symmetric(D_wb[1:2*interior_indmax,1:2*interior_indmax])

    display("Diagonalizing D")
    eigenmodes = diagonalize_D(D)

    create_group(analysis_file, "eigenmodes")
    analysis_file["eigenmodes"]["eigvals"] = eigenmodes["eigvals"]

    #Keep first 16 modes
    analysis_file["eigenmodes"]["eigvecs"] = size(eigenmodes["eigvecs"])[2]<=16 ? eigenmodes["eigvecs"] :  eigenmodes["eigvecs"][:,1:16]
    x0int = extract_frame_data_for_type("x",1,frames_support[string(length(frames_support))])
    y0int = extract_frame_data_for_type("y",1,frames_support[string(length(frames_support))])
    analysis_file["eigenmodes"]["x0int"] = x0int #for visualization purposes
    analysis_file["eigenmodes"]["y0int"] = y0int #for visualization purposes

    create_group(analysis_file, "projs")
    v_projs = project_on_eigvecs(eigenmodes["eigvecs"], vx,vy)
    analysis_file["projs"]["v_projs"] = v_projs



    FT_v_projs = temporal_Fourier_transform(dt, v_projs, min_t_ind = min_t_ind, output_not_avg=true)
    save_dict2h5!(analysis_file,FT_v_projs, "FT_v_projs")




    return analysis_file
end




function readdir_filt(folder_path)
    directories = filter!(e->e!=".DS_Store",readdir(folder_path))
    return directories

end

function save_dict2h5!(file,dict, dictsavename)

    create_group(file, dictsavename)

    for (key, val) in dict
        file[dictsavename][string(key)] = val
    end
    return file
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



function extract_frame_data_for_type(datakey, type, frame_data)

    return frame_data[datakey][ frame_data["type"].==type ]

end

