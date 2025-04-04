include(joinpath("..","src","Engine.jl"))
include("AnalysisPipeline.jl")


base_folder = "/data1/kammeraat/sa/phi_1/Nlin_20/vary_J_Dr/" 
analysis_base_folder = joinpath(base_folder, "analysis")

plot_base_folder = mkpath(joinpath(base_folder, "plots")) 

tree = construct_folder_tree_param_param_seed(analysis_base_folder)


using PDFmerger
using CairoMakie

#Individual plots
function plot_mean_phi_over_time(seed,seedanalysis_file)
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

    append_pdf!( joinpath(plot_base_folder,"mean_phi_over_time.pdf"), "temp.pdf", cleanup=true)
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

    append_pdf!( joinpath(plot_base_folder,"K_over_time.pdf"), "temp.pdf", cleanup=true)
end
acces_param1_param2_seedanalysis(tree, [plot_K_over_time])
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

    append_pdf!( joinpath(plot_base_folder,"mean_px_time.pdf"), "temp.pdf", cleanup=true)
end

acces_param1_param2_seedanalysis(tree, [plot_mean_phi_over_time,plot_px_over_time,plot_psi_over_time])

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
    xlims!(0,Δt[1:length(Cavg)][end])
    lines!(ax, Δt[1:length(Cavg)],Cavg)
    display(f)

    save("temp.pdf",f)

    append_pdf!( joinpath(plot_base_folder,"AUTO_p.pdf"), "temp.pdf", cleanup=true)
end
acces_param1_param2_seedanalysis(tree, [plot_AUTO_p])

function plot_v_p_projections(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    v_projs = seedanalysis_file["v_projs"]
    p_projs = seedanalysis_file["p_projs"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    ax = Axis(f[1,1], xlabel="v_proj", ylabel=" p_proj", title="J=$J, Dr=$Dr")
    scatter!(ax, v_projs[1,:],p_projs[1,:], color=:blue)

    scatter!(ax, v_projs[2,:],p_projs[2,:], color=:red)

    display(f)

    save("temp.pdf",f)

    append_pdf!( joinpath(plot_base_folder,"v_p_projs.pdf"), "temp.pdf", cleanup=true)
end
acces_param1_param2_seedanalysis(tree, [plot_v_p_projections])

function plot_p_projections(seed,seedanalysis_file)
    f = Figure()
    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    v0 = seedanalysis_file["v0"]
    p_projs = seedanalysis_file["p_projs"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    ax = Axis(f[1,1], xlabel="mode number", ylabel="proj^2", title="J=$J, Dr=$Dr", yscale=log10)
    scatter!(ax, mean(p_projs.^2, dims=2)[:,1], color=:blue)

    display(f)

    save("temp.pdf",f)

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
    ax = Axis(f[1,1], xlabel="mode number", ylabel="dis_proj^2", title="J=$J, Dr=$Dr", yscale=log10)
    scatter!(ax, mean(dis_projs.^2, dims=2)[:,1], color=:blue)

    display(f)

    save("temp.pdf",f)

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
    ax = Axis(f[1,1], xlabel="ω", ylabel="dis_proj^2", title="J=$J, Dr=$Dr", yscale=log10)
    scatter!(ax,ωs, mean(dis_projs[:,500:end].^2, dims=2)[:,1], color=:blue)

    display(f)

    save("temp.pdf",f)

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
    ax = Axis(f[1,1], xlabel="ω", ylabel="v_proj^2", title="J=$J, Dr=$Dr")#, yscale=log10)
    scatter!(ax,ωs, mean(v_projs[:,500:end].^2, dims=2)[:,1], color=:blue)

    display(f)

    save("temp.pdf",f)

    append_pdf!( joinpath(plot_base_folder,"lin_omega_v_projs.pdf"), "temp.pdf", cleanup=true)
end


acces_param1_param2_seedanalysis(tree, [plot_ω_v_projections])
acces_param1_param2_seedanalysis(tree, [plot_ω_dis_projections])
acces_param1_param2_seedanalysis(tree, [plot_dis_projections])
acces_param1_param2_seedanalysis(tree, [plot_p_projections])
#### Combined plots
Drs = []
Js = []
ϕms = []
ψms = []
ψmstds = []

Kms = []
Kmstds = []
function collect_params(seed,seedanalysis_file)

    Dr = seedanalysis_file["Dr"]
    J = seedanalysis_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
    t = seedanalysis_file["integration_info"]["save_tax"]
    push!(Drs, Dr)
    push!(Js, J)

    push!(ϕms, mean(seedanalysis_file["mean_ϕ"][t .> 500]))

    push!(ψms, mean(seedanalysis_file["ψ"][t .> 500]))

    push!(ψmstds, std(seedanalysis_file["ψ"][t .> 500]))

    push!(Kms, mean(seedanalysis_file["K"][t .> 500]))

    push!(Kmstds, std(seedanalysis_file["K"][t .> 500]))

end

acces_param1_param2_seedanalysis(tree, [collect_params])

f = Figure();
ax = Axis(f[1,1], xlabel="Dr", ylabel="J", xscale=log10);
hm=heatmap!(Drs, Js, ϕms)

Colorbar(f[:, end+1], hm,label=" ⟨ ⟨ϕ⟩ ⟩")

scatter!(ax, Drs, Js, color=:white, strokecolor=:black, strokewidth=1)

display(f);
save(joinpath(plot_base_folder,"mean2_phi_Dr_J.pdf"),f)

f = Figure();
ax = Axis(f[1,1], xlabel="Dr", ylabel="J", xscale=log10);
hm=heatmap!(Drs, Js, abs.(ϕms))

Colorbar(f[:, end+1], hm,label=" |⟨ ⟨ϕ⟩ ⟩|")

scatter!(ax, Drs, Js, color=:white, strokecolor=:black, strokewidth=1)

display(f);
save(joinpath(plot_base_folder,"abs_mean2_phi_Dr_J.pdf"),f)


f = Figure();
ax = Axis(f[1,1], xlabel="Dr", ylabel="J", xscale=log10);
hm=heatmap!(Drs, Js, ψms)

Colorbar(f[:, end+1], hm,label=" ⟨ ψ(t) ≡ 1/N |∑ e^(iϕ)| ⟩_t ")

scatter!(ax, Drs, Js, color=:white, strokecolor=:black, strokewidth=1)

display(f);
save(joinpath(plot_base_folder,"mean_psi_Dr_J.pdf"),f)


f = Figure();
ax = Axis(f[1,1], xlabel="Dr", ylabel="J", xscale=log10);
hm=heatmap!(Drs, Js, Kms)

Colorbar(f[:, end+1], hm,label=" ⟨ K(t) ≡ 1/N |∑ e^(iθ_p)| ⟩_t ")

scatter!(ax, Drs, Js, color=:white, strokecolor=:black, strokewidth=1)

display(f);
save(joinpath(plot_base_folder,"mean_K_Dr_J.pdf"),f)



f = Figure();
ax = Axis(f[1,1], xlabel="1/Dr", ylabel="⟨ ψ(t)⟩_t", xscale=log10);
ylims!(ax, (0,1))
colormap=:viridis
for i in eachindex(Drs)

    Dr = Drs[i]
    J = Js[i]
    ψm = ψms[i]
    ψmstd = ψmstds[i]

    if J==1. || J==0.1 || true
        if Dr==0.01
            scatter!(ax, 1/Dr, ψm, color=J, colorrange=(0,1.), label="$J", colormap=colormap)
        else
            scatter!(ax, 1/Dr, ψm, color=J, colorrange=(0,1.),colormap=colormap)
        end
        errorbars!(ax,[1/Dr],[ψm], [ψmstd], color=[J], colorrange=(0,1.), colormap=colormap)
    end

end
f[1,2]=Legend(f,ax, "J")
display(f)
