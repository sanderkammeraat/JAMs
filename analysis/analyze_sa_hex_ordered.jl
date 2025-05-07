

include(joinpath("..","src","Engine.jl"))
include("AnalysisPipeline.jl")
include("AnalysisFunctions.jl")

#base_folder = "/data1/kammeraat/sa/phi_1/Nlin_20/vary_J_Dr"

base_folder = joinpath(homedir(),"sa","survey","hex_ordered","phi_1","Nlin_4","vary_J_Dr")

raw_data_base_folder = joinpath(base_folder, "simdata")

analysis_base_folder = mkpath(joinpath(base_folder, "analysis_FT"))

#Make tree to navigate simulation data folder structure
tree = construct_folder_tree_param_param_seed(raw_data_base_folder)

function analyze_single_seed_inner!(analysis_file, system, integration_info, frames; frames_support) 

    t = integration_info["save_tax"]

    #For reference


    Nt = length(t)
    #Ignore boundary particles
    Np = length(extract_frame_data_for_type("id",1,frames["1"]))

    # collect x, first index is fastest

    x = zeros(Np, Nt)
    y = zeros(Np, Nt)

    px = zeros(Np, Nt)
    py = zeros(Np, Nt)
    vx = zeros(Np, Nt)
    vy = zeros(Np, Nt)
    θp = zeros(Np, Nt)
    θv = zeros(Np, Nt)
    ϕ = zeros(Np, Nt)
    ψ = zeros(Nt)
    #Kuramoto phase coherence
    K = zeros(Nt)

    @views for i in 1:Nt
        x[:,i] = extract_frame_data_for_type("x", 1, frames[string(i)])
        y[:,i] = extract_frame_data_for_type("y", 1, frames[string(i)])



        vx[:,i] = extract_frame_data_for_type("vx", 1, frames[string(i)])
        vy[:,i] = extract_frame_data_for_type("vy", 1, frames[string(i)])

        px[:,i] = extract_frame_data_for_type("px", 1, frames[string(i)])
        py[:,i] = extract_frame_data_for_type("py", 1, frames[string(i)])

        θv[:,i] = angle.(vx[:,i] .+ 1im .* vy[:,i]) 
        θp[:,i] = angle.(px[:,i] .+ 1im .* py[:,i]) 
        phase_factor = exp.(1im .* (θv[:,i] - θp[:,i]))
        ϕ[:,i] = angle.( phase_factor)

        ψ[i] = 1/Np*norm(sum(phase_factor))

        K[i] = 1/Np*norm(sum(exp.(1im .* θp[:,i])))
    end

    #Set maximum value of radial bin edge
    xmax = maximum(x)
    ymax = maximum(y)

    rmax = 2*sqrt(xmax^2 + ymax^2)



    #Include boundary points!
    x0 = frames["1"]["x"]
    y0 = frames["1"]["y"]


    
    x0int = x[:,1]
    y0int = y[:,1]




    k = system["forces"]["pair_forces"]["soft_disk_force"]["karray"]
    R = frames["1"]["R"]
    type = frames["1"]["type"]


    dis_x = zeros(Np, Nt)
    dis_y = zeros(Np, Nt)

    @views for i in 1:Nt
        dis_x[:,i] = x[:,i] .- x[:,1]
        dis_y[:,i] = y[:,i] .- y[:,1]

    end



    #With boundary
    D_wb = construct_D(x0, y0, k, R, type)

    #Boundary particles are stored at the end
    interior_indmax = sum(type.==1)
    D = Symmetric(D_wb[1:2*interior_indmax,1:2*interior_indmax])


    

    modes = diagonalize_D(D)

    #save_dict!(analysis_file, modes, "modes")

    #dis_projs = project_on_eigvecs(modes["eigvecs"], dis_x,dis_y)

    #p_projs = project_on_eigvecs(modes["eigvecs"], px,py)

    

    #analysis_file["p_projs"] = p_projs



    AUTO_p = auto_correlation(t, px, py, minrow=500)


    save_dict!(analysis_file,secondary_temporal_Fourier_transform(AUTO_p["Δt"][2]-AUTO_p["Δt"][1], AUTO_p["Cavg"]), "FT_auto_P")


    #save_dict!(analysis_file,spatiotemporal_p_correlation(2.5, rmax, x[:,1],y[:,1],px, py, min_t_ind=500), "SPTE_p" )
    save_dict!(analysis_file, AUTO_p, "AUTO_p")

    dt = t[2] - t[1]

    save_dict!(analysis_file,temporal_Fourier_transform(dt, dis_x;  min_t_ind=500), "FT_dx")

    save_dict!(analysis_file,temporal_Fourier_transform(dt, px;  min_t_ind=500), "FT_px")

    save_dict!(analysis_file,temporal_Fourier_transform(dt, vx;  min_t_ind=500), "FT_vx")



    save_dict!(analysis_file, spatial_p_correlation(2.5, rmax, x,y,px, py,min_t_ind=500), "SPAT_p" )

    v_projs = project_on_eigvecs(modes["eigvecs"], vx,vy)

    #analysis_file["dis_projs"] = dis_projs

    analysis_file["v_projs"] = v_projs

    analysis_file["R"] = frames["1"]["R"]
    analysis_file["id"] = frames["1"]["id"]
    analysis_file["type"] = frames["1"]["type"]
    analysis_file["v0"] = frames["1"]["v0"][1]
    analysis_file["Dr"] = frames["1"]["Dr"][1]

    #analysis_file["D"] = D

    analysis_file["modes/eigvals"] = modes["eigvals"]


    #only store the first, at most 100, eigvecs
    analysis_file["modes/100eigvecs"] =  size(modes["eigvecs"])[2]<=100 ? modes["eigvecs"] :  modes["eigvecs"][:,1:100]

    analysis_file["x0"] = x0
    analysis_file["y0"] = y0
    
    analysis_file["x0int"] = x0int
    analysis_file["y0int"] = y0int

    #analysis_file["θv"] = θv   
    #analysis_file["θp"] = θp  
    #analysis_file["ϕ"] = ϕ

    #analysis_file["ψ"] = ψ

    analysis_file["mean_ψ"] = mean(ψ, dims=1)[1,:] 

    analysis_file["std_ψ"] = std(ψ, dims=1)[1,:] 

    #analysis_file["K"] = K

    analysis_file["mean_K"] = mean(K, dims=1)[1,:] 

    analysis_file["std_K"] = std(K, dims=1)[1,:] 



    analysis_file["mean_ϕ"] =mean(ϕ, dims=1)[1,:] 

    analysis_file["std_ϕ"] =std(ϕ, dims=1)[1,:] 

    analysis_file["mean_px"] =mean(px, dims=1)[1,:] 

    analysis_file["std_px"] =std(px, dims=1)[1,:] 

    analysis_file["t"] = integration_info["save_tax"]


    return analysis_file

end

run_multithreaded_analysis_param1_param2_seed(tree, analyze_single_seed_inner!, analysis_base_folder, overwrite=true)
