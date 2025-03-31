include(joinpath("..","src","Engine.jl"))
include("AnalysisPipeline.jl")



analysis_base_folder = joinpath(homedir(),"sa", "varyJ", "analysis")

plot_base_folder = mkpath(joinpath(homedir(),"sa", "varyJ", "plots")) 

tree = construct_folder_tree_param_param_seed(analysis_base_folder)


using PDFmerger
using CairoMakie


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

acces_param1_param2_seedanalysis(tree, [plot_mean_phi_over_time])



