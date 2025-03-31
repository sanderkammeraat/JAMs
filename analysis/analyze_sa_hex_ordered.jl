

include(joinpath("..","src","Engine.jl"))
include("AnalysisPipeline.jl")
include("AnalysisFunctions.jl")

raw_data_base_folder = joinpath(homedir(),"sa", "varyJ", "simdata")

analysis_base_folder = mkpath(joinpath(homedir(),"sa", "varyJ", "analysis"))

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

    @views for i in 1:Nt
        for j in 1:Np
            x[:,i] = extract_frame_data_for_type("x", 1, frames[string(i)])
            y[:,i] = extract_frame_data_for_type("y", 1, frames[string(i)])

            vx[:,i] = extract_frame_data_for_type("vx", 1, frames[string(i)])
            vy[:,i] = extract_frame_data_for_type("vy", 1, frames[string(i)])

            px[:,i] = extract_frame_data_for_type("px", 1, frames[string(i)])
            py[:,i] = extract_frame_data_for_type("py", 1, frames[string(i)])

            θv[:,i] = angle.(vx[:,i] .+ 1im .* vy[:,i]) 
            θp[:,i] = angle.(px[:,i] .+ 1im .* py[:,i]) 
            ϕ[:,i] = angle.( exp.(1im .* (θv[:,i] - θp[:,i])))
        end
    end

    analysis_file["x"] = x
    analysis_file["y"] = y

    analysis_file["θv"] = θv   
    analysis_file["θp"] = θp  
    analysis_file["ϕ"] = ϕ


    analysis_file["mean_ϕ"] =mean(ϕ, dims=1)[1,:] 

    analysis_file["std_ϕ"] =std(ϕ, dims=1)[1,:] 


    return analysis_file

end

run_serial_analysis_param1_param2_seed(tree, analyze_single_seed_inner!, analysis_base_folder, overwrite=true)

#original = jldopen( joinpath(raw_data_base_folder,"Dr_0.01","J_0.0","seed_1",raw_data_file_name ), "r")
#test = jldopen( joinpath(analysis_base_folder,"Dr_0.01","J_0.0","seed_1.jld2"), "r")
#close(test)


##Develop room for inner function

#Check with single 

raw_data_file_path =  joinpath(raw_data_base_folder,"Dr_0.01","J_1.0","seed_8","raw_data.jld2")


analysis_file_path = joinpath( mkpath(joinpath(analysis_base_folder,"Dr_0.01","J_1.0")), "seed_8.jld2")


analyze_single_seed_outer(raw_data_file_path, analysis_file_path, analyze_single_seed_inner!, overwrite=true)


result = jldopen(analysis_file_path,"r")
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


close(result)