

include(joinpath("..","src","Engine.jl"))
include("AnalysisPipeline.jl")
include("AnalysisFunctions.jl")


#base_folder = joinpath(homedir(),"sa","survey","hex_disordered","phi_1","Nlin_4","vary_J_Dr")

#base_folder = "/data1/kammeraat/sa/survey/hex_disordered/phi_1/Nlin_20/vary_J_Dr/" 

base_folder = joinpath(homedir(),"sa","survey","hex_disordered","phi_1","Nlin_4","vary_J_Dr")
raw_data_base_folder = joinpath(base_folder, "simdata")

analysis_base_folder = mkpath(joinpath(base_folder, "analysis_FT"))

#Make tree to navigate simulation data folder structure
tree = construct_folder_tree_param_param_seed(raw_data_base_folder)

function analyze_single_seed_inner!(analysis_file, system, integration_info, frames; frames_support) 

    t = integration_info["save_tax"]



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



    #Include boundary points!
    x0 = frames_support[string(length(frames_support))]["x"]
    y0 = frames_support[string(length(frames_support))]["y"]

    x0int = extract_frame_data_for_type("x",1,frames_support[string(length(frames_support))])
    y0int = extract_frame_data_for_type("y",1,frames_support[string(length(frames_support))])



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
    #analysis_file["D"] = D

    analysis_file["modes/eigvals"] = modes["eigvals"]


    #only store the first, at most 100, eigvecs
    analysis_file["modes/100eigvecs"] =  size(modes["eigvecs"])[2]<=100 ? modes["eigvecs"] :  modes["eigvecs"][:,1:100]

    #dis_projs = project_on_eigvecs(modes["eigvecs"], dis_x,dis_y)

    #p_projs = project_on_eigvecs(modes["eigvecs"], px,py)

    v_projs = project_on_eigvecs(modes["eigvecs"], vx,vy)

    #analysis_file["dis_projs"] = dis_projs

    analysis_file["v_projs"] = v_projs
    
    #analysis_file["v_proj_2_tmean"] = mean(v_projs[:,500:end].^2, dims=2)[:,1]

    #analysis_file["p_projs"] = p_projs

    #Set maximum value of radial bin edge
    xmax = maximum(x)
    ymax = maximum(y)

    rmax = 2*sqrt(xmax^2 + ymax^2)

    #save_dict!(analysis_file,spatiotemporal_p_correlation(2.5, rmax, x0int,y0int,px, py, min_t_ind=500), "SPTE_p" )


    AUTO_p = auto_correlation(t, px, py, minrow=500)


    save_dict!(analysis_file,secondary_temporal_Fourier_transform(AUTO_p["Δt"][2]-AUTO_p["Δt"][1], AUTO_p["Cavg"]), "FT_AUTO_p")


    #save_dict!(analysis_file,spatiotemporal_p_correlation(2.5, rmax, x[:,1],y[:,1],px, py, min_t_ind=500), "SPTE_p" )
    save_dict!(analysis_file, AUTO_p, "AUTO_p")

    dt = t[2] - t[1]

    save_dict!(analysis_file,temporal_Fourier_transform(dt, dis_x;  min_t_ind=500), "FT_dx")

    save_dict!(analysis_file,temporal_Fourier_transform(dt, px;  min_t_ind=500), "FT_px")

    save_dict!(analysis_file,temporal_Fourier_transform(dt, vx;  min_t_ind=500), "FT_vx")



    

    save_dict!(analysis_file, spatial_p_correlation(2.5, rmax, x,y,px, py,min_t_ind=500), "SPAT_p" )


    #For reference
    analysis_file["R"] = frames["1"]["R"]
    analysis_file["id"] = frames["1"]["id"]
    analysis_file["type"] = frames["1"]["type"]
    analysis_file["v0"] = frames["1"]["v0"][1]
    analysis_file["Dr"] = frames["1"]["Dr"][1]


    analysis_file["x0"] = x0
    analysis_file["y0"] = y0

    analysis_file["x0int"] = x0int
    analysis_file["y0int"] = y0int



    #analysis_file["θv"] = θv   
    #analysis_file["θp"] = θp  
    #analysis_file["ϕ"] = ϕ

    #analysis_file["ψ"] = ψ

    analysis_file["mean_ψ"] = mean(ψ)

    analysis_file["std_ψ"] = std(ψ)

    #analysis_file["K"] = K

    analysis_file["mean_K"] = mean(K)

    analysis_file["std_K"] = std(K)



    analysis_file["mean_ϕ"] =mean(ϕ, dims=1)[1,:] 

    analysis_file["std_ϕ"] =std(ϕ, dims=1)[1,:] 

    analysis_file["mean_px"] =mean(px, dims=1)[1,:] 

    analysis_file["std_px"] =std(px, dims=1)[1,:] 

    analysis_file["t"] = integration_info["save_tax"]


    return analysis_file

end



for (key,val) in tree
    println(key)
end

partial_tree=Dict()


partial_tree["J_0.0"] = tree["J_0.0"]

partial_tree["J_0.5"] = tree["J_0.5"]

partial_tree["J_5.0"] = tree["J_5.0"]


run_multithreaded_analysis_param1_param2_seed(tree, analyze_single_seed_inner!, analysis_base_folder, overwrite=true, raw_data_file_name="sa_raw_data.jld2", support_raw_data_file_name="ra_raw_data.jld2")

