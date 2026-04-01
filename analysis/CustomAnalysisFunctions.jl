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

    J = haskey(system["forces"]["external"],"self_align_with_v_unit_force") ? system["forces"]["external"]["self_align_with_v_unit_force"]["β"] : system["forces"]["external"]["self_align_with_v_force"]["β"]
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
    #
    v_projs = project_on_eigvecs(eigenmodes["eigvecs"], vx,vy)
    analysis_file["projs"]["v_projs"] = v_projs

    FT_v_projs = temporal_Fourier_transform(dt, v_projs, min_t_ind = min_t_ind, output_not_avg=true)
    save_dict2h5!(analysis_file, FT_v_projs, "FT_v_projs")
    #

    #
    p_projs = project_on_eigvecs(eigenmodes["eigvecs"], px,py)
    analysis_file["projs"]["p_projs"] = p_projs



    FT_p_projs = temporal_Fourier_transform(dt, p_projs, min_t_ind = min_t_ind, output_not_avg=true)
    save_dict2h5!(analysis_file, FT_p_projs, "FT_p_projs")
    #


    run_sa_analysis_add_auto_p!(analysis_file, raw_data_file, support_raw_data_file = support_raw_data_file)
    run_sa_analysis_add_auto_v!(analysis_file, raw_data_file, support_raw_data_file = support_raw_data_file)

    run_sa_analysis_add_spatial_cor!(analysis_file, raw_data_file, support_raw_data_file = support_raw_data_file)



    
    return analysis_file
end

function run_sa_analysis_add_auto_p!(analysis_file, raw_data_file; support_raw_data_file = nothing)

    if !haskey(analysis_file,"auto_p")
        frames = raw_data_file["frames"]

        system = raw_data_file["system"]

        frames_support = support_raw_data_file["frames"]

        t =  raw_data_file["integration_info"]["save_tax"]

        v0 = frames["1"]["v0"][1]

        Dr = frames["1"]["Dr"][1]

        J = system["forces"]["external"]["self_align_with_v_unit_force"]["β"]

        k = system["forces"]["pair"]["soft_disk_force"]["karray"]

        R = frames["1"]["R"]

        type = frames["1"]["type"]

        Nt = length(t)

        #Ignore boundary particles
        Nint = length(extract_frame_data_for_type("id",1,frames["1"]))


        min_t_ind = 500


        dt = t[2] - t[1]

        px = zeros(Nint, Nt)
        py = zeros(Nint, Nt)


        @views for i in 1:Nt

            px[:,i] .= extract_frame_data_for_type("px", 1, frames[string(i)])
            py[:,i] .= extract_frame_data_for_type("py", 1, frames[string(i)])

        end

        auto_p = auto_correlation(t[end-1000:end],px[:,end-1000:end], py[:,end-1000:end], minrow=1)

        save_dict2h5!(analysis_file, auto_p, "auto_p")
    end
    return analysis_file
    GC.gc()
end

function run_sa_analysis_add_auto_v!(analysis_file, raw_data_file; support_raw_data_file = nothing)

    if !haskey(analysis_file,"auto_v")
        frames = raw_data_file["frames"]

        system = raw_data_file["system"]

        frames_support = support_raw_data_file["frames"]

        t =  raw_data_file["integration_info"]["save_tax"]

        v0 = frames["1"]["v0"][1]

        Dr = frames["1"]["Dr"][1]

        J = system["forces"]["external"]["self_align_with_v_unit_force"]["β"]

        k = system["forces"]["pair"]["soft_disk_force"]["karray"]

        R = frames["1"]["R"]

        type = frames["1"]["type"]

        Nt = length(t)

        #Ignore boundary particles
        Nint = length(extract_frame_data_for_type("id",1,frames["1"]))


        


        dt = t[2] - t[1]

        vx = zeros(Nint, Nt)
        vy = zeros(Nint, Nt)


        @views for i in 1:Nt

            vx[:,i] .= extract_frame_data_for_type("vx", 1, frames[string(i)])
            vy[:,i] .= extract_frame_data_for_type("vy", 1, frames[string(i)])

        end


        auto_v = auto_correlation(t[end-1000:end],vx[:,end-1000:end], vy[:,end-1000:end], minrow=1,normalized=true)

        save_dict2h5!(analysis_file, auto_v, "auto_v")
    end
    return analysis_file
    GC.gc()
end

function run_sa_analysis_add_spatial_cor!(analysis_file, raw_data_file; support_raw_data_file = nothing)

    if !haskey(analysis_file,"spatial_vcor") || !haskey(analysis_file,"spatial_pcor") 
        frames = raw_data_file["frames"]

        system = raw_data_file["system"]

        frames_support = support_raw_data_file["frames"]

        t =  raw_data_file["integration_info"]["save_tax"]

        v0 = frames["1"]["v0"][1]

        Dr = frames["1"]["Dr"][1]

        J = system["forces"]["external"]["self_align_with_v_unit_force"]["β"]

        k = system["forces"]["pair"]["soft_disk_force"]["karray"]

        R = frames["1"]["R"]

        type = frames["1"]["type"]

        Nt = length(t)

        #Ignore boundary particles
        Nint = length(extract_frame_data_for_type("id",1,frames["1"]))

        min_t_ind = 500

        dt = t[2] - t[1]

        x = zeros(Nint, Nt)
        y = zeros(Nint, Nt)

        vx = zeros(Nint, Nt)
        vy = zeros(Nint, Nt)

        px = zeros(Nint, Nt)
        py = zeros(Nint, Nt)

        @views for i in 1:Nt

            x[:,i] .= extract_frame_data_for_type("x", 1, frames[string(i)])
            y[:,i] .= extract_frame_data_for_type("y", 1, frames[string(i)])

            px[:,i] .= extract_frame_data_for_type("px", 1, frames[string(i)])
            py[:,i] .= extract_frame_data_for_type("py", 1, frames[string(i)])

            vx[:,i] .= extract_frame_data_for_type("vx", 1, frames[string(i)])
            vy[:,i] .= extract_frame_data_for_type("vy", 1, frames[string(i)])

        end

        x0 = frames_support[string(length(frames_support))]["x"]
        y0 = frames_support[string(length(frames_support))]["y"]

        nc = 10

        spatial_vcor = spatial_v_correlation(3, 200, x[:,min_t_ind:nc:end], y[:,min_t_ind:nc:end], vx[:,min_t_ind:nc:end], vy[:,min_t_ind:nc:end])

        save_dict2h5!(analysis_file, spatial_vcor, "spatial_vcor")

        spatial_pcor = spatial_p_correlation(3, 200, x[:,min_t_ind:nc:end], y[:,min_t_ind:nc:end], px[:,min_t_ind:nc:end], py[:,min_t_ind:nc:end])

        save_dict2h5!(analysis_file, spatial_pcor, "spatial_pcor")


    end
    return analysis_file
    GC.gc()
end



function run_free_sa_analysis!(analysis_file, raw_data_file; support_raw_data_file = nothing)


    frames = raw_data_file["frames"]

    system = raw_data_file["system"]

    frames_support = support_raw_data_file["frames"]

    t =  raw_data_file["integration_info"]["save_tax"]
    analysis_file["t"] = t

    v0 = frames["1"]["v0"][1]
    analysis_file["v0"] = v0

    Dr = frames["1"]["Dr"][1]
    analysis_file["Dr"] = Dr

    J = haskey(system["forces"]["external"],"self_align_with_v_unit_force") ? system["forces"]["external"]["self_align_with_v_unit_force"]["β"] : system["forces"]["external"]["self_align_with_v_force"]["β"]
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
    D_wb = construct_D(x0, y0, k, R, type, periodic_system_sizes = system["sizes"][1:2])


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

    v_projs = project_on_eigvecs(eigenmodes["eigvecs"], vx , vy)
    analysis_file["projs"]["v_projs"] = v_projs


    FT_v_projs = temporal_Fourier_transform(dt, v_projs, min_t_ind = min_t_ind, output_not_avg=true)
    save_dict2h5!(analysis_file, FT_v_projs, "FT_v_projs")


    run_sa_analysis_add_auto_p!(analysis_file, raw_data_file, support_raw_data_file = support_raw_data_file)


    return analysis_file
end
function run_free_BP_analysis!(analysis_file, raw_data_file; support_raw_data_file = nothing)


    frames = raw_data_file["frames"]

    system = raw_data_file["system"]

    frames_support = support_raw_data_file["frames"]

    t =  raw_data_file["integration_info"]["save_tax"]
    analysis_file["t"] = t

    v0 = frames["1"]["v0"][1]
    analysis_file["v0"] = v0

    Dr = frames["1"]["Dr"][1]
    analysis_file["Dr"] = Dr


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



    #Heavy work: doing the mode analysis
    #With boundary
    display("Constructing D_wb")
    D = construct_D(x0, y0, k, R, type, periodic_system_sizes = system["sizes"][1:3])

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

    v_projs = project_on_eigvecs(eigenmodes["eigvecs"], vx , vy)
    analysis_file["projs"]["v_projs"] = v_projs


    return analysis_file
end



function correct_vrms_analysis!(analysis_file, raw_data_file; support_raw_data_file = nothing)


    frames = raw_data_file["frames"]

    system = raw_data_file["system"]

    frames_support = support_raw_data_file["frames"]

    t =  raw_data_file["integration_info"]["save_tax"]
    analysis_file["t"] = t

    v0 = frames["1"]["v0"][1]
    analysis_file["v0"] = v0

    Dr = frames["1"]["Dr"][1]
    analysis_file["Dr"] = Dr

    J = haskey(system["forces"]["external"],"self_align_with_v_unit_force") ? system["forces"]["external"]["self_align_with_v_unit_force"]["β"] : system["forces"]["external"]["self_align_with_v_force"]["β"]
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


    @views for i in 1:Nt
        x[:,i] .= extract_frame_data_for_type("x", 1, frames[string(i)])
        y[:,i] .= extract_frame_data_for_type("y", 1, frames[string(i)])

        vx[:,i] .= extract_frame_data_for_type("vx", 1, frames[string(i)])
        vy[:,i] .= extract_frame_data_for_type("vy", 1, frames[string(i)])

    end

    #vrms info 

    vrms_particle_time_avg = sqrt.(mean((vx.^2 + vy.^2)[:,min_t_ind:end],dims=2)[:,1])
    analysis_file["vrms_particle_time_avg"] = vrms_particle_time_avg

    vrms_particle_avg_time_avg = sqrt.( mean( (vx.^2 + vy.^2)[:,min_t_ind:end] ))
    analysis_file["vrms_particle_avg_time_avg"] =  vrms_particle_avg_time_avg


    vrms_particle_avg_time = sqrt.( mean((vx.^2 + vy.^2)[:,min_t_ind:end],dims=1)[1,:] )
    analysis_file["vrms_particle_avg_time"] = vrms_particle_avg_time


    
    return analysis_file
end

function run_cs_analysis!(analysis_file, raw_data_file; support_raw_data_file = nothing)

    frames = raw_data_file["frames"]
    system = raw_data_file["system"]
    t =  raw_data_file["integration_info"]["save_tax"]
    analysis_file["t"] = t

    v0 = frames["1"]["v0"][1]
    analysis_file["v0"] = v0

    Dr = frames["1"]["Dr"][1]
    analysis_file["Dr"] = Dr


    k = system["forces"]["pair"]["soft_disk_force"]["karray"]
    analysis_file["k"] = k

    analysis_file["morse_Dearray"] = system["forces"]["pair"]["morse_force"]["Dearray"]
    analysis_file["morse_aarray"] = system["forces"]["pair"]["morse_force"]["aarray"]

    R = frames["1"]["R"]
    analysis_file["R"] = R


    type = frames["1"]["type"]
    analysis_file["type"] = type


    Nt = length(t)
    analysis_file["Nt"] = Nt
    #Ignore boundary particles
    Nact = length(extract_frame_data_for_type("id",2,frames["1"]))
    analysis_file["Nact"] = Nact

    Npas = length(extract_frame_data_for_type("id",1,frames["1"]))
    analysis_file["Npas"] = Npas

    min_t_ind = 1
    analysis_file["min_t_ind"] = min_t_ind

    dt = t[2] - t[1]
    analysis_file["dt"] = dt


    x = zeros(Nact, Nt)
    y = zeros(Nact, Nt)

    xuw = zeros(Nact, Nt)
    yuw = zeros(Nact, Nt)

    vx = zeros(Nact, Nt)
    vy = zeros(Nact, Nt)

    px = zeros(Nact, Nt)
    py = zeros(Nact, Nt)

    #qx = zeros(Nint, Nt)
    #qy = zeros(Nint, Nt)


    @views for i in 1:Nt
        x[:,i] .= extract_frame_data_for_type("x", 2, frames[string(i)])
        y[:,i] .= extract_frame_data_for_type("y", 2, frames[string(i)])

        xuw[:,i] .= extract_frame_data_for_type("xuw", 2, frames[string(i)])
        yuw[:,i] .= extract_frame_data_for_type("yuw", 2, frames[string(i)])

        vx[:,i] .= extract_frame_data_for_type("vx", 2, frames[string(i)])
        vy[:,i] .= extract_frame_data_for_type("vy", 2, frames[string(i)])

        px[:,i] .= extract_frame_data_for_type("px", 2, frames[string(i)])
        py[:,i] .= extract_frame_data_for_type("py", 2, frames[string(i)])

        #qx[:,i] .= extract_frame_data_for_type("qx", 1, frames[string(i)])
        #qy[:,i] .= extract_frame_data_for_type("qy", 1, frames[string(i)])

    end

    vrms_particle_time_avg = sqrt.(mean((vx.^2 + vy.^2)[:,min_t_ind:end],dims=2)[:,1])
    analysis_file["vrms_particle_time_avg"] = vrms_particle_time_avg

    vrms_particle_avg_time_avg = sqrt.( mean( (vx.^2 + vy.^2)[:,min_t_ind:end] ))
    analysis_file["vrms_particle_avg_time_avg"] =  vrms_particle_avg_time_avg


    vrms_particle_avg_time = sqrt.( mean((vx.^2 + vy.^2)[:,min_t_ind:end],dims=1)[1,:] )
    analysis_file["vrms_particle_avg_time"] = vrms_particle_avg_time

    msd = MSD(t,xuw,yuw)
    save_dict2h5!(analysis_file, msd, "MSD")

    auto_p = auto_correlation(t[min_t_ind:end],px[:,min_t_ind:end], py[:,min_t_ind:end], minrow=1)
    save_dict2h5!(analysis_file, auto_p, "auto_p")

    auto_v_norm = auto_correlation(t[min_t_ind:end],vx[:,min_t_ind:end], vy[:,min_t_ind:end], minrow=1,normalized=true)

    save_dict2h5!(analysis_file, auto_v_norm, "auto_v_norm")

    return analysis_file
end

