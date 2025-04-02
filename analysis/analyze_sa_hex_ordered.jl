

include(joinpath("..","src","Engine.jl"))
include("AnalysisPipeline.jl")
include("AnalysisFunctions.jl")

raw_data_base_folder = joinpath(homedir(),"sa", "vary_J_Dr", "simdata")

analysis_base_folder = mkpath(joinpath(homedir(),"sa", "vary_J_Dr", "analysis"))

#Make tree to navigate simulation data folder structure
tree = construct_folder_tree_param_param_seed(raw_data_base_folder)

function analyze_single_seed_inner!(analysis_file, system, integration_info, frames) 

    t = integration_info["save_tax"]

    #For reference
    analysis_file["R"] = frames["1"]["R"]
    analysis_file["id"] = frames["1"]["id"]
    analysis_file["type"] = frames["1"]["type"]
    analysis_file["v0"] = frames["1"]["v0"][1]
    analysis_file["Dr"] = frames["1"]["Dr"][1]

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

    analysis_file["x"] = x
    analysis_file["y"] = y

    analysis_file["AUTO_p"] = auto_correlation(t, px, py, minrow=500)

    #Include boundary points!
    x0 = frames["1"]["x"]
    y0 = frames["1"]["y"]
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


    analysis_file["D"] = D

    modes = diagonalize_D(D)
    analysis_file["modes"] = modes

    dis_projs = project_on_eigvecs(modes["eigvecs"], dis_x,dis_y)

    p_projs = project_on_eigvecs(modes["eigvecs"], px,py)

    v_projs = project_on_eigvecs(modes["eigvecs"], vx,vy)

    analysis_file["dis_projs"] = dis_projs

    analysis_file["v_projs"] = v_projs

    analysis_file["p_projs"] = p_projs

    analysis_file["θv"] = θv   
    analysis_file["θp"] = θp  
    analysis_file["ϕ"] = ϕ

    analysis_file["ψ"] = ψ

    analysis_file["mean_ψ"] = mean(ψ)

    analysis_file["std_ψ"] = std(ψ)

    analysis_file["K"] = K

    analysis_file["mean_K"] = mean(K)

    analysis_file["std_K"] = std(K)



    analysis_file["mean_ϕ"] =mean(ϕ, dims=1)[1,:] 

    analysis_file["std_ϕ"] =std(ϕ, dims=1)[1,:] 

    analysis_file["mean_px"] =mean(px, dims=1)[1,:] 

    analysis_file["std_px"] =std(px, dims=1)[1,:] 


    return analysis_file

end

run_serial_analysis_param1_param2_seed(tree, analyze_single_seed_inner!, analysis_base_folder, overwrite=true)

#original = jldopen( joinpath(raw_data_base_folder,"Dr_0.01","J_0.0","seed_1",raw_data_file_name ), "r")
#test = jldopen( joinpath(analysis_base_folder,"Dr_0.01","J_0.0","seed_1.jld2"), "r")
#close(test)


##Develop room for inner function

#Check with single 


raw_data_file_path =  joinpath(raw_data_base_folder,"Dr_0.01","J_1.0","seed_50","raw_data.jld2")


analysis_file_path = joinpath( mkpath(joinpath(analysis_base_folder,"Dr_0.01","J_1.0")), "seed_50.jld2")


analyze_single_seed_outer(raw_data_file_path, analysis_file_path, analyze_single_seed_inner!, overwrite=true)


result = jldopen(analysis_file_path,"r")


plot(result["type"])
display(result)

θv = result["θv"]
θp = result["θp"]
ϕ = result["ϕ"]
f = Figure();
ax = Axis(f[1,1])
for i in 500:600

    scatter!(ax, θv[:,i], θp[:,i])
end

display(f)

ax2 = Axis(f[1,2])

lines!(ax2, result["integration_info"]["save_tax"],mean(ϕ, dims=1)[1,:])

display(f)


AUTO_p = result["AUTO_p"]

Δt = AUTO_p["Δt"]


heatmap(AUTO_p["C"])



lines(AUTO_p["Cavg"])

close(result)


