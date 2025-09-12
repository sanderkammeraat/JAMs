include(joinpath("..","src","Engine.jl"))

include("AnalysisPipeline.jl")


base_folder = joinpath("/Volumes","T7_Shield","sa","survey","hex_disordered", "phi_1", "Nlin_20", "vary_J_Dr")

begin

#base_folder = joinpath("/data1/kammeraat/sa/phi_1/Nlin_20/vary_J_Dr/" 

#base_folder = joinpath(homedir(),"mounting","data1_kammeraat","sa", "phi_1", "Nlin_20", "vary_J_Dr")

#base_folder = joinpath(homedir(), "sa", "survey","hex_disordered","phi_1", "Nlin_4", "vary_J_Dr")



analysis_base_folder = joinpath(base_folder, "analysis_FT")

plot_base_folder = mkpath(joinpath(base_folder, "plots_11_9")) 

tree = construct_folder_tree_param_param_seed(analysis_base_folder)


using PDFmerger
using CairoMakie

end 

function get_tag(seedanalysis_file)

    v0 = v0 = seedanalysis_file["v0"]

    type = seedanalysis_file["type"]

    R = seedanalysis_file["R"]

    k = seedanalysis_file["system"]["forces"]["pair_forces"]["soft_disk_force"]["karray"]

    #Interior particles
    Nint = sum(type .== 1)
    ϕ = 1
    #Check if all radii are the same, if so do, else , hardcoded 0.15, because I did not store the polydispersity of the initial conditions in the  analysis file
    poly =  all( R .== R[1]) ? 0. : 0.15


    tag = Dict("ϕ"=>ϕ, "v0"=> v0, "Nint"=> Nint, "poly"=>poly, "k"=>k)

    return tag

end



#Individual plots
function plot_phi_over_time(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    ϕm  = seedanalysis_file["mean_ϕ"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    ax = Axis(f[1,1], xlabel="t", ylabel=" ⟨ϕ⟩", title="J=$J, Dr=$Dr")
    lines!(ax, t, ϕm)
    display(f)

    save("temp.pdf",f)
    
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"phi_over_time.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"phi_over_time.pdf"), "temp.pdf", cleanup=true)
end

function plot_psi_over_time(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    ψ  = seedanalysis_file["ψ"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    ax = Axis(f[1,1], xlabel="t", ylabel=" ψ(t)", title="J=$J, Dr=$Dr")
    lines!(ax, t, ψ)
    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"psi_over_time.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"psi_over_time.pdf"), "temp.pdf", cleanup=true)
end

function plot_K_over_time(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    K  = seedanalysis_file["K"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    ax = Axis(f[1,1], xlabel="t", ylabel=" K(t)", title="J=$J, Dr=$Dr")
    ylims!(0,1)
    lines!(ax, t, K)
    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"K_over_time.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"K_over_time.pdf"), "temp.pdf", cleanup=true)
end

function plot_px_over_time(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    pxm  = seedanalysis_file["mean_px"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    ax = Axis(f[1,1], xlabel="t", ylabel=" ⟨p_x⟩", title="J=$J, Dr=$Dr")
    lines!(ax, t,pxm)
    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"px_time.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"px_time.pdf"), "temp.pdf", cleanup=true)
end




function plot_AUTO_p(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    AUTO_p = seedanalysis_file["AUTO_p"]
    Cavg  = AUTO_p["Cavg"]
    Δt  = AUTO_p["Δt"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    ax = Axis(f[1,1], xlabel="Δt", ylabel=" ⟨ p(t+Δt) ⋅ p(t) ⟩_{t>500}", title="J=$J, Dr=$Dr")
    ylims!(-1,1)
    xlims!(0,500)#Δt[1:length(Cavg)][end])
    lines!(ax, Δt[1:length(Cavg)],Cavg)
    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"AUTO_p.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"AUTO_p.pdf"), "temp.pdf", cleanup=true)
end


function plot_v_p_projections(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    v_projs = seedanalysis_file["v_projs"]
    p_projs = seedanalysis_file["p_projs"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    ax = Axis(f[1,1], xlabel="⟨λ_n|V⟩ ", ylabel="⟨λ_n|P⟩", title="J=$J, Dr=$Dr")
    scatter!(ax, v_projs[1,:],p_projs[1,:], color=:blue)

    scatter!(ax, v_projs[2,:],p_projs[2,:], color=:red)

    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"v_p_projs.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"v_p_projs.pdf"), "temp.pdf", cleanup=true)
end


function plot_p_projections(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    p_projs = seedanalysis_file["p_projs"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    ax = Axis(f[1,1], xlabel="n", ylabel="⟨λ_n|P⟩^2", title="J=$J, Dr=$Dr", yscale=log10)
    scatter!(ax, mean(p_projs.^2, dims=2)[:,1], color=:blue)

    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"p_projs.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"p_projs.pdf"), "temp.pdf", cleanup=true)
end

function plot_dis_projections(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    dis_projs = seedanalysis_file["dis_projs"]
    ωs = sqrt.(seedanalysis_file["modes"]["eigvals"])
    t = seedanalysis_file["integration_info"]["save_tax"]
    ax = Axis(f[1,1], xlabel="n", ylabel="⟨λ_n|δR⟩^2", title="J=$J, Dr=$Dr", yscale=log10)
    scatter!(ax, mean(dis_projs.^2, dims=2)[:,1], color=:blue)

    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"dis_projs.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"dis_projs.pdf"), "temp.pdf", cleanup=true)
end

function plot_ω_dis_projections(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    dis_projs = seedanalysis_file["dis_projs"]
    ωs = sqrt.(seedanalysis_file["modes"]["eigvals"])

    t = seedanalysis_file["integration_info"]["save_tax"]
    ax = Axis(f[1,1], xlabel="ω_n", ylabel="⟨λ_n=ω_n^2|δR⟩^2", title="J=$J, Dr=$Dr", yscale=log10)
    scatter!(ax,ωs, mean(dis_projs[:,500:end].^2, dims=2)[:,1], color=:blue)

    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"omega_dis_projs.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"omega_dis_projs.pdf"), "temp.pdf", cleanup=true)
end

function plot_ω_v_projections(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    v_projs = seedanalysis_file["v_projs"]
    ωs = sqrt.(seedanalysis_file["modes"]["eigvals"])
    
    t = seedanalysis_file["integration_info"]["save_tax"]
    
    ax = Axis(f[1,1], xlabel="ω_n", ylabel="⟨λ_n=ω_n^2|V⟩^2", title="J=$J, Dr=$Dr")#, yscale=log10)
    scatter!(ax,ωs, mean(v_projs[:,500:end].^2, dims=2)[:,1], color=:blue)
    #ylims!(1e-7, 1e-5)
    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"omega_v_projs.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"omega_v_projs.pdf"), "temp.pdf", cleanup=true)
end


function plot_t_v_projections(seed,seedanalysis_file)

    

    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    #if J==0.1
    with_theme(theme_latexfonts()) do 
    f = Figure()


    
    
    v_projs = seedanalysis_file["v_projs"]
    #
    
    t = seedanalysis_file["integration_info"]["save_tax"]
    
    ax = Axis(f[1,1], xlabel="t", ylabel=L"Vel. proj.: $\langle \lambda_n|\delta \dot{R} \rangle$", title="J=$J, Dr=$Dr")#, yscale=log10)

    marker_labels=[
        (:circle, ":circle"),
        (:rect, ":rect"),
        (:diamond, ":diamond"),
        (:hexagon, ":hexagon"),
        (:cross, ":cross"),
        (:xcross, ":xcross"),
        (:utriangle, ":utriangle"),
        (:dtriangle, ":dtriangle"),
        (:ltriangle, ":ltriangle"),
        (:rtriangle, ":rtriangle"),
        (:pentagon, ":pentagon")]


    Nmodes = 5

    for i in 1:Nmodes
        scatterlines!(ax,t, v_projs[i,:], color=i, colorrange=(1, Nmodes), label="$i", marker=marker_labels[i][1])

    end

    tag = get_tag(seedanalysis_file)

    # lines!(ax,t, v_projs, color=:blue)

    f[1,2]=Legend(f,ax, L"Mode numbers")
    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
    xlims!(1000, 5000)
    #ylims!(-0.03, 0.03)
    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"t_v_projs.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"t_v_projs.pdf"), "temp.pdf", cleanup=true)
end
end
begin
    collective_plot_file_name ="t_v_projs.pdf"
    try 
        rm(joinpath(plot_base_folder,collective_plot_file_name))
    catch
    end
    acces_param1_param2_seedanalysis(tree, [plot_t_v_projections])
    end

function plot_FT_px(seed,seedanalysis_file)

    

    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    #if J==0.1
    with_theme(theme_latexfonts()) do 
    f = Figure()

    eigvals = seedanalysis_file["modes"]["eigvals"]
    
    N = length(eigvals)/2

    theory = sqrt( J*v0/(2*N) * sum(    1 ./ ( 1 ./ eigvals .+ 1/Dr)   ))


    v_projs = seedanalysis_file["v_projs"]

    theory_empirical =  sqrt( J/(N*v0) * sum(  eigvals .* mean(v_projs[:,500:end].^2, dims=2)[:,1]) / sqrt(sum(mean(v_projs[:,500:end].^2, dims=2)[:,1]))  )

    
    
    FT = seedanalysis_file["FT_px"]

    
    t = seedanalysis_file["integration_info"]["save_tax"]


    ax = Axis(f[1,1], yscale=log10, xlabel=L"\omega", ylabel=L"|  \mathcal{F}\{p_x(t)\}(\omega)|^2", title="J=$J, Dr=$Dr", xscale=log10);
    scatter!(ax, FT["ω"], FT["pavg_X2"])
    scatter!(ax,FT["ω_max"] , FT["max_X2"],label="max", color="orange")

    #vlines!(ax,theory,label="theory 1st order J", color="green")

    #vlines!(ax,theory_empirical,label="theory empirical projection", color="pink")

    tag = get_tag(seedanalysis_file)

    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

    f[1,2]=Legend(f,ax)
    ylims!(ax, low=1e1,high= 1e7)

    xlims!(ax, low=1e-4,high= 3)
    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"FT_px.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"FT_px.pdf"), "temp.pdf", cleanup=true)
end
end
begin
    collective_plot_file_name ="FT_px.pdf"
    try 
        rm(joinpath(plot_base_folder,collective_plot_file_name))
    catch
    end
    acces_param1_param2_seedanalysis(tree, [plot_FT_px])
end

function plot_FT_vx(seed,seedanalysis_file)
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    #if J==0.1
    with_theme(theme_latexfonts()) do 
    f = Figure()

    FT = seedanalysis_file["FT_vx"]
    
    t = seedanalysis_file["integration_info"]["save_tax"]

    ax = Axis(f[1,1], yscale=log10, xlabel=L"\omega", ylabel=L"|  \mathcal{F}\{v_x(t)\}(\omega)|^2", title="J=$J, Dr=$Dr");
    scatter!(ax, FT["ω"], FT["pavg_X2"])
    scatter!(ax,FT["ω_max"] , FT["max_X2"],label="max")

    tag = get_tag(seedanalysis_file)

    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

    f[1,2]=Legend(f,ax)
    #ylims!(ax, low=1e-7,high= 1e2)
    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"FT_vx.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"FT_vx.pdf"), "temp.pdf", cleanup=true)
end
end
begin
    collective_plot_file_name ="FT_vx.pdf"
    try 
        rm(joinpath(plot_base_folder,collective_plot_file_name))
    catch
    end
    acces_param1_param2_seedanalysis(tree, [plot_FT_vx])
end
function plot_FT_dx(seed,seedanalysis_file)
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    #if J==0.1
    with_theme(theme_latexfonts()) do 
    f = Figure()

    FT = seedanalysis_file["FT_dx"]
    
    t = seedanalysis_file["integration_info"]["save_tax"]

    ax = Axis(f[1,1], yscale=log10, xlabel=L"\omega", ylabel=L"|  \mathcal{F}\{δx(t)\}(\omega)|^2", title="J=$J, Dr=$Dr");
    scatter!(ax, FT["ω"], FT["pavg_X2"])
    scatter!(ax,FT["ω_max"] , FT["max_X2"],label="max")

    tag = get_tag(seedanalysis_file)

    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

    f[1,2]=Legend(f,ax)
    #ylims!(ax, low=1e-7,high= 1e2)
    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"FT_dx.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"FT_dx.pdf"), "temp.pdf", cleanup=true)
end
end
begin
    collective_plot_file_name ="FT_dx.pdf"
    try 
        rm(joinpath(plot_base_folder,collective_plot_file_name))
    catch
    end
    acces_param1_param2_seedanalysis(tree, [plot_FT_dx])
end
function plot_FT_AUTO_p(seed,seedanalysis_file)
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    #if J==0.1
    with_theme(theme_latexfonts()) do 
    f = Figure()

    FT = seedanalysis_file["FT_auto_P"]
    
    t = seedanalysis_file["integration_info"]["save_tax"]

    ax = Axis(f[1,1], yscale=log10, xlabel=L"\omega", ylabel=L"|  \mathcal{F}\{AUTO_p(Δt)\}(\omega)|^2", title="J=$J, Dr=$Dr");
    scatter!(ax, FT["ω"], FT["X2"])
    scatter!(ax,FT["ω_max"] , FT["max_X2"],label="max")

    tag = get_tag(seedanalysis_file)

    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

    f[1,2]=Legend(f,ax)
    #ylims!(ax, low=1e-7,high= 1e2)
    display(f)

    save("temp.pdf",f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"FT_AUTO_p.pdf"),f)

    append_pdf!( joinpath(plot_base_folder,"FT_AUTO_p.pdf"), "temp.pdf", cleanup=true)
end
end
begin
    collective_plot_file_name ="FT_AUTO_p.pdf"
    try 
        rm(joinpath(plot_base_folder,collective_plot_file_name))
    catch
    end
    acces_param1_param2_seedanalysis(tree, [plot_FT_AUTO_p])
end

#acces_param1_param2_seedanalysis(tree, [plot_v_p_projections])


#acces_param1_param2_seedanalysis(tree, [plot_ω_v_projections])
#acces_param1_param2_seedanalysis(tree, [plot_ω_dis_projections])
#acces_param1_param2_seedanalysis(tree, [plot_dis_projections])
#acces_param1_param2_seedanalysis(tree, [plot_p_projections])
#### Combined plots
begin
Drs = []
Js = []
ϕms = []

pxms = []

ψms = []
ψmstds = []

Kms = []
Kmstds = []

ωmaxs = []

tags=[]
function collect_params(seed,seedanalysis_file)

    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    if Dr>0 && J>0
        push!(Drs, Dr)
        push!(Js, J)
        t_ind_transient = 500

        push!(ϕms, mean((seedanalysis_file["mean_ϕ"])[t_ind_transient:end]))

        push!(pxms, mean((seedanalysis_file["mean_px"])[t_ind_transient:end]))
        FT = seedanalysis_file["FT_px"]
        push!(ωmaxs, FT["ω_max"])

        #push!(ψms, (seedanalysis_file["mean_ψ"])[t_ind_transient:end])

        #push!(ψmstds, seedanalysis_file["std_ψ"][t_ind_transient:end])

        #println(seedanalysis_file["mean_K"])
        #push!(Kms, mean( ( seedanalysis_file["mean_K"])[t_ind_transient:end] ) )

        #push!(Kmstds, seedanalysis_file["std_K"][t_ind_transient:end])

        push!(tags,get_tag(seedanalysis_file))
    end
end

acces_param1_param2_seedanalysis(tree, [collect_params])

function Dr_J_heatmap(Drs, Js, vals, cbarlabel,tags,  save_path)


    f = Figure();
    ax = Axis(f[1,1], xlabel="log_10(Dr)", ylabel="log_10(J)");


    hm=heatmap!(ax,log10.(Drs), log10.(Js), vals)

    Colorbar(f[:, end+1], hm,label=cbarlabel)

    scatter!(ax, log10.(Drs), log10.(Js), color=:white, strokecolor=:black, strokewidth=1)


    tag = tags[1]
    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

    display(f);

    save(save_path,f)
end

save_path = joinpath(plot_base_folder,"phi_Dr_J.pdf")
Dr_J_heatmap(Drs, Js, ϕms, "⟨ ⟨ϕ⟩ ⟩",tags, save_path)



save_path = joinpath(plot_base_folder,"abs_phi_Dr_J.pdf")
Dr_J_heatmap(Drs, Js, abs.(ϕms), "|⟨ ⟨ϕ⟩ ⟩|",tags, save_path)


save_path = joinpath(plot_base_folder,"px_omega_max_Dr_J.pdf")
Dr_J_heatmap(Drs, Js, log10.(ωmaxs), "log_10(ω_max(px))",tags, save_path)
end

#save_path = joinpath(plot_base_folder,"psi_Dr_J.pdf")
#Dr_J_heatmap(Drs, Js, ψms, "⟨ ψ(t) ≡ 1/N |∑ e^(iϕ)| ⟩_t ",tags, save_path)

#save_path = joinpath(plot_base_folder,"K_Dr_J.pdf")
#Dr_J_heatmap(Drs, Js, Kms, " ⟨ K(t) ≡ 1/N |∑ e^(iθ_p)| ⟩_t ",tags, save_path)



# f = Figure();
# ax = Axis(f[1,1], xlabel="1/Dr", ylabel="⟨ ψ(t)⟩_t", xscale=log10);
# ylims!(ax, (0,1))
# colormap=:viridis
# for i in eachindex(Drs)

#     Dr = Drs[i]
#     J = Js[i]
#     ψm = ψms[i]
#     ψmstd = ψmstds[i]

#     if J==1. || J==0.1 || true
#         if Dr==0.01
#             scatter!(ax, 1/Dr, ψm, color=J, colorrange=(0,1.), label="$J", colormap=colormap)
#         else
#             scatter!(ax, 1/Dr, ψm, color=J, colorrange=(0,1.),colormap=colormap)
#         end
#         errorbars!(ax,[1/Dr],[ψm], [ψmstd], color=[J], colorrange=(0,1.), colormap=colormap)
#     end

# end
# f[1,2]=Legend(f,ax, "J")
# display(f)

#J with different Drs
begin
    collective_plot_file_name ="J_omega_v_proj.pdf"
    try 
        rm(joinpath(plot_base_folder,collective_plot_file_name))
    catch
    end
    for (param1,_) in sort(tree)


        with_theme(theme_latexfonts()) do 
            f = Figure()
            ax = Axis(f[1,1], xlabel=L"ω_n", ylabel= L"Vel. proj.: $\langle \lambda_n|\delta \dot{R} \rangle^2$", title="$param1", yscale=log10)#, xscale=log10)

            
            marker_ind=1

            marker_labels=[
                (:circle, ":circle"),
                (:rect, ":rect"),
                (:diamond, ":diamond"),
                (:hexagon, ":hexagon"),
                (:cross, ":cross"),
                (:xcross, ":xcross"),
                (:utriangle, ":utriangle"),
                (:dtriangle, ":dtriangle"),
                (:ltriangle, ":ltriangle"),
                (:rtriangle, ":rtriangle"),
                (:pentagon, ":pentagon")]
            for (param2,_) in sort(tree[param1])
                for (seed, seedpath) in tree[param1][param2]
                    jldopen(seedpath, "r") do seedanalysis_file
                        Dr = seedanalysis_file["Dr"]
                        
                        J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
                        v0 = seedanalysis_file["v0"]
                        v_projs = seedanalysis_file["v_projs"]
                        if all(seedanalysis_file["modes"]["eigvals"].>0)
                        ωs = sqrt.(seedanalysis_file["modes"]["eigvals"])
                        
                        t = seedanalysis_file["integration_info"]["save_tax"]
                        if Dr!=0
                            tau =1/Dr
                            theory = v0^2  ./ (2 .+ 2 .* ωs.^2 .* tau)
                            lines!(ax, ωs, theory ,color=marker_ind, colorrange=(1,length(tree[param1])),colormap=Reverse(:gist_rainbow), alpha=1, linestyle=:dash)


                            
                        end
                        numerics =  mean(v_projs[:,500:end].^2, dims=2)[:,1]
                        scatterlines!(ax,ωs, numerics,  label="$Dr",colormap=Reverse(:gist_rainbow),color=marker_ind, colorrange=(1,length(tree[param1])), linewidth=0.5, marker=marker_labels[marker_ind][1], )
                        marker_ind+=1
                        ylims!(ax,low=5e-9,high=1e-1)
                        xlims!(ax,low=0,high=2.5)
                        global tag = get_tag(seedanalysis_file)
                        end
                    end
                end
            end
        f[1,2]=Legend(f,ax, L"D_r")
        Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
        display(f)
        save("temp.pdf",f)
        append_pdf!( joinpath(plot_base_folder,collective_plot_file_name), "temp.pdf", cleanup=true)
        subfolder_path = mkpath(joinpath(plot_base_folder, "$param1"))
        save(joinpath(subfolder_path,"omega_v_proj.pdf"),f )
        end
    end
end

begin
    collective_plot_file_name ="J_n_v_proj.pdf"
    try 
        rm(joinpath(plot_base_folder,collective_plot_file_name))
    catch
    end
    for (param1,_) in sort(tree)


        with_theme(theme_latexfonts()) do 
            f = Figure()
            ax = Axis(f[1,1], xlabel=L"n", ylabel= L"Vel. proj.: $\langle \lambda_n|\delta \dot{R} \rangle^2$", title="$param1", yscale=log10)#, xscale=log10)

            
            marker_ind=1

            marker_labels=[
                (:circle, ":circle"),
                (:rect, ":rect"),
                (:diamond, ":diamond"),
                (:hexagon, ":hexagon"),
                (:cross, ":cross"),
                (:xcross, ":xcross"),
                (:utriangle, ":utriangle"),
                (:dtriangle, ":dtriangle"),
                (:ltriangle, ":ltriangle"),
                (:rtriangle, ":rtriangle"),
                (:pentagon, ":pentagon")]
            for (param2,_) in sort(tree[param1])
                for (seed, seedpath) in tree[param1][param2]
                    jldopen(seedpath, "r") do seedanalysis_file
                        Dr = seedanalysis_file["Dr"]
                        
                        J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
                        v0 = seedanalysis_file["v0"]
                        v_projs = seedanalysis_file["v_projs"]
                        if all(seedanalysis_file["modes"]["eigvals"].>0)
                        ωs = sqrt.(seedanalysis_file["modes"]["eigvals"])
                        
                        t = seedanalysis_file["integration_info"]["save_tax"]
                        numerics =  mean(v_projs[:,500:end].^2, dims=2)[:,1]
                        scatterlines!(ax,1:length(ωs), numerics,  label="$Dr",colormap=Reverse(:gist_rainbow),color=marker_ind, colorrange=(1,length(tree[param1])), marker=marker_labels[marker_ind][1], linewidth=0.1)
                        if Dr!=0
                            tau =1/Dr
                            theory = v0^2  ./ (2 .+ 2 .* ωs.^2 .* tau)
                            lines!(ax, 1:length(ωs), theory ,color=marker_ind, colorrange=(1,length(tree[param1])),colormap=Reverse(:gist_rainbow), alpha=1)
                        end
                        
                        marker_ind+=1
                        ylims!(ax,low=5e-9,high=1e-1)
                        xlims!(ax,low=0,high=40)
                        global tag = get_tag(seedanalysis_file)
                        end
                    end
                end
            end
        f[1,2]=Legend(f,ax, L"D_r")
        Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
        display(f)
        save("temp.pdf",f)
        append_pdf!( joinpath(plot_base_folder,collective_plot_file_name), "temp.pdf", cleanup=true)
        subfolder_path = mkpath(joinpath(plot_base_folder, "$param1"))
        save(joinpath(subfolder_path,"n_v_proj.pdf"),f )
        end
    end
end



begin
    collective_plot_file_name ="J_auto_p.pdf"
    try 
        rm(joinpath(plot_base_folder,collective_plot_file_name))
    catch
    end
    for (param1,_) in sort(tree)
        with_theme(theme_latexfonts()) do 
        f = Figure()
        ax = Axis(f[1,1], xlabel=L"t", ylabel= L"$\langle \mathbf{p}(t_0+t) \cdot \mathbf{p}(t_0) \rangle$", title="$param1")#, xscale=log10)
    
        
        marker_ind=1
    
        marker_labels=[
            (:circle, ":circle"),
            (:rect, ":rect"),
            (:diamond, ":diamond"),
            (:hexagon, ":hexagon"),
            (:cross, ":cross"),
            (:xcross, ":xcross"),
            (:utriangle, ":utriangle"),
            (:dtriangle, ":dtriangle"),
            (:ltriangle, ":ltriangle"),
            (:rtriangle, ":rtriangle"),
            (:pentagon, ":pentagon")]
        for (param2,_) in sort(tree[param1])
            for (seed, seedpath) in tree[param1][param2]
                jldopen(seedpath, "r") do seedanalysis_file
    
    
                    Dr = seedanalysis_file["Dr"]
                    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
                    v0 = seedanalysis_file["v0"]
                    AUTO_p = seedanalysis_file["AUTO_p"]
                    
                    Δt  =AUTO_p["Δt"]
                    C = AUTO_p["Cavg"]
                    
                    scatterlines!(ax,Δt[1:length(C)], C,  label="$Dr",colormap=Reverse(:gist_rainbow),color=marker_ind, colorrange=(1,length(tree[param1])), linewidth=0.5, marker=marker_labels[marker_ind][1], )
                    marker_ind+=1
                    ylims!(ax, low=-1.01,high=1.01)
                    xlims!(ax, low=0, high=100)
    
                    global tag = get_tag(seedanalysis_file)
                end
            end
        end
    
        f[1,2]=Legend(f,ax, L"D_r")
        Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
        display(f)
        save("temp.pdf",f)
        append_pdf!( joinpath(plot_base_folder,collective_plot_file_name), "temp.pdf", cleanup=true)
        subfolder_path = mkpath(joinpath(plot_base_folder, "$param1"))
        save(joinpath(subfolder_path,"auto_p.pdf"),f )
    end
    end
end





begin
    collective_plot_file_name ="J_spatial_p.pdf"
    try 
        rm(joinpath(plot_base_folder,collective_plot_file_name))
    catch
    end
    for (param1,_) in sort(tree)
    
    
        with_theme(theme_latexfonts()) do 
        f = Figure()
        ax = Axis(f[1,1], xlabel=L"r", ylabel= L"$\langle \mathbf{p}(t_0, \mathbf{r_0}) \cdot \mathbf{p}(t_0, \mathbf{r_0}+\mathbf{r})) \rangle$", title="$param1")#, xscale=log10)
    
        
        marker_ind=1
    
        marker_labels=[
            (:circle, ":circle"),
            (:rect, ":rect"),
            (:diamond, ":diamond"),
            (:hexagon, ":hexagon"),
            (:cross, ":cross"),
            (:xcross, ":xcross"),
            (:utriangle, ":utriangle"),
            (:dtriangle, ":dtriangle"),
            (:ltriangle, ":ltriangle"),
            (:rtriangle, ":rtriangle"),
            (:pentagon, ":pentagon")]
        for (param2,_) in sort(tree[param1])
            for (seed, seedpath) in tree[param1][param2]
                jldopen(seedpath, "r") do seedanalysis_file
    
    
                    Dr = seedanalysis_file["Dr"]
                    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
                    v0 = seedanalysis_file["v0"]
                    SPAT_p = seedanalysis_file["SPAT_p"]
                    
                    rbc  =SPAT_p["rbc"]
                    C = SPAT_p["C"]


                    Cavg = zeros(size(C)[2])
                    for i in 1:size(C)[2]
                        Cavg[i] = mean( filter!( e->!isnan(e) ,C[:,i] ) )

                    end
                    
                    scatterlines!(ax,rbc, Cavg,  label="$Dr",colormap=Reverse(:gist_rainbow),color=marker_ind, colorrange=(1,length(tree[param1])), linewidth=0.5, marker=marker_labels[marker_ind][1], )
                    marker_ind+=1
                    ylims!(ax, low=-1.01,high=1.01)
                    #xlims!(ax, low=0, high=100)
    
                    global tag = get_tag(seedanalysis_file)
                end
            end
    
        end
        f[1,2]=Legend(f,ax, L"D_r")
        Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
        display(f)
        save("temp.pdf",f)
        append_pdf!( joinpath(plot_base_folder,collective_plot_file_name), "temp.pdf", cleanup=true)
        subfolder_path = mkpath(joinpath(plot_base_folder, "$param1"))
        save(joinpath(subfolder_path,"spatial_p.pdf"),f )
    end
    end 
end




function plot_first_n_modes(seed,seedanalysis_file)

    with_theme(theme_latexfonts()) do 
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    pxm  = seedanalysis_file["mean_px"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    f = Figure(size=(2000,2000));
    mm=1

    x0  = seedanalysis_file["x0"]
    y0  = seedanalysis_file["y0"]

    R = seedanalysis_file["R"]

    id = seedanalysis_file["id"]
    
    n = 3
    for i in 1:n
        for j in 1:n
        print(mm)
    
        m = mm

        eigvecs = seedanalysis_file["modes"]["100eigvecs"]
        eigvals = seedanalysis_file["modes"]["eigvals"]
    
        eigen_m_value = eigvals[m]
        eigen_m_vector =  reshape(eigvecs[:,m], (2, round(Int64,length(eigvecs[:,m])/2)) )
        eigen_m_x = eigen_m_vector[1,:]
        eigen_m_y = eigen_m_vector[2,:]
    
    
        ax=Axis(f[i,j], xlabel = "x", ylabel="y", aspect=DataAspect(), title=L"λ_%$(m)= %$(eigen_m_value)")
        scatter!(ax,x0,y0, markerspace=:data, markersize=2 .*R, color = id,marker = Circle,alpha=0.3)
    
        scale = length(R)>100 ? 100 : 4


        arrows!(ax, x0, y0, eigen_m_x, eigen_m_y, lengthscale = scale)
        mm+=1
        end
    
    end

    tag = get_tag(seedanalysis_file)
    Label(f[n+1,1],"J = $J, Dr=$(Dr), "*"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save(joinpath(subfolder_path,"mode_vis.pdf"),f)

    display(f)
end

end

acces_param1_param2_seedanalysis(tree, [plot_first_n_modes])


function plot_spatiotemporal_p(seed,seedanalysis_file)
    with_theme(theme_latexfonts()) do 

    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    SPTE_p = seedanalysis_file["SPTE_p"]
    ax = Axis(f[1,1], yreversed=true, xlabel=L"r", ylabel=L"t")
    t_max_ind = 5000

    t = seedanalysis_file["t"]
    
    hm=heatmap!(ax,  SPTE_p["rbe"],t[1:t_max_ind],   transpose(SPTE_p["C"][1:t_max_ind,:]), colormap=:seismic, colorrange=(-1,1))
    vlines!(ax, SPTE_p["rbe"], color=:black)

    tag = get_tag(seedanalysis_file)
    Label(f[end+1,1],"J = $J, Dr=$(Dr), "*"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
    Colorbar(f[1, end+1], hm,label=L"$\langle \mathbf{p}(t_s+t, \mathbf{r_0}) \cdot \mathbf{p}(t_s, \mathbf{r_0}+\mathbf{r})) \rangle$")
    display(f)
    subfolder_path = mkpath(joinpath(plot_base_folder, "J_$J", "Dr_$Dr"))

    save("temp.pdf",f)
    append_pdf!( joinpath(plot_base_folder,"spatiotemporal.pdf"), "temp.pdf", cleanup=true)


    
    save(joinpath(subfolder_path,"spatiotemporal_p.pdf"),f)
end


end




acces_param1_param2_seedanalysis(tree, [plot_spatiotemporal_p])

