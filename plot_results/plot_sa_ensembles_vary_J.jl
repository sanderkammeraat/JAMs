
using PDFmerger
include("../analysis/AnalysisPipeline.jl")

#base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/Nlin_20"

base_folder = "/Users/kammeraat/mounting/data2_kammeraat/sa/statistics/hex_disordered/phi_1.3/Nlin_20"

figure_save_folder = joinpath(base_folder, "figures_24_11")
mkpath(figure_save_folder)



ef  = findfile(base_folder, "ensemble.h5")


ensemble_files = [ load_file(ef[i]) for i in eachindex(ef)]

Drs = [e["Dr"] for e in ensemble_files]


Js = [e["J"] for e in ensemble_files]

v0 = [e["v0"] for e in ensemble_files] #by color

GLMakie.activate!()
begin

f = Figure()

ax = Axis(f[1,1], xlabel=L"r", ylabel=L"vrms(r)/v0")

for e in ensemble_files
    Je = e["J"]

    scatterlines!(ax, e["vrms_r"]["r_bin_centers"], e["vrms_r"]["vrms_r"]/e["v0"], colorrange = (0, 1 ) , color=e["J"], label="J = $(e["J"])")

end
ylims!(ax,0,1.5)
f[1,2]=Legend(f,ax)
display(f)

end

begin


f = Figure()

ax = Axis(f[1,1], xlabel=L"λ", ylabel=L"v projs /v_0^2", yscale=log10)

for e in ensemble_files


    J= e["J"]
    Dr = e["Dr"]
    v0 = e["v0"]

    

        display("plottings")
        tau =1/Dr

        eigval_bin_centers = e["v_projs_time_avg"]["eigval_bin_centers"]

        eigval_ind = 1 

        theory_ABP = v0^2  ./ (2 .+ 2 .* eigval_bin_centers .* tau)  /v0^2


        a_min = sqrt(1 +  (eigval_bin_centers[eigval_ind]/(2 * J) + 1/(2 * tau *J ))^2 ) - (eigval_bin_centers[eigval_ind]/(2 * J) + 1/(2 * tau *J))
        a=a_min
        #a based on loweest mode selection
        #a = e["vrms"]/v0
        #display(a)
        a = e["vrms"]/v0

        the_eigvals = collect(range(0, stop=5, step=0.001))

        A =  sqrt.( (1/tau + J * a) .* the_eigvals)

        B = the_eigvals .+  .-J /a  .+ 1/tau  .+ J*a

        select = B .>= 2* A

        A =   A[select]
        B = B[select]
        theory_amin = @. pi * ( B - sqrt(B^2 - 4 * A^2))/B/sqrt(-4*A^2 + 2 * B * (B - sqrt(B^2-4*A^2)) ) * 2/tau /4/pi 

        scatterlines!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["J"], colorrange = (0, 1) ,  label="J = $(e["J"])")

        lines!(ax,eigval_bin_centers,theory_ABP, color=e["J"], colorrange = (0, 1) ,  label="J = $(e["J"])")
        #lines!(ax,the_eigvals[select],theory_amin, color=e["J"], colorrange = (0, 1) ,  label="J = $(e["J"])")

 
end
xlims!(ax, 0,1)
f[1,2]=Legend(f,ax)
display(f)
end







begin #heatmap wn vs w
using CairoMakie
CairoMakie.activate!()
for e in ensemble_files

    f = Figure()

    J= e["J"]
    Dr = e["Dr"]
    v0 = e["v0"]

    ax = Axis(f[1,1], ylabel=L"ω_n", xlabel=L"ω", title="J=$J, Dr = $Dr, v0 = $(v0)")

    tau =1/Dr

    eigval_bin_centers = e["FT_v_projs"]["eigval_bin_centers"]

    w = e["FT_v_projs"]["w"]

    X = e["FT_v_projs"]["X2"]

    heatmap!(ax,w, sqrt.(eigval_bin_centers), log10.(transpose(X)), rasterize=true, colorrange=(-5,5))
    cb = Colorbar(f[1, 2], limits = (-5, 5), colormap = :viridis, flipaxis = false, label = "log10(ω^2 |a_n(ω)|^2)")


    xlims!(0,1)
    ylims!(0,2.5)
    # theory_ABP = v0^2  ./ (2 .+ 2 .* eigval_bin_centers .* tau)  /v0^2

    # scatterlines!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["v0"], colorrange = (0, 0.01 ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")

    # lines!(ax,eigval_bin_centers,theory_ABP, colorrange = (0, maximum(Drs) ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")
    #f[1,2]=Legend(f,ax)
    save("temp.pdf",f)
    append_pdf!( joinpath(figure_save_folder,"wn_w.pdf"), "temp.pdf", cleanup=true)

    display(f)
end
end
begin #heatmap wn vs w THEORY
using CairoMakie
CairoMakie.activate!()
for e in ensemble_files

    f = Figure()

    J= e["J"]
    Dr = e["Dr"]
    v0 = e["v0"]

    vrms = e["vrms"]
    #a = vrms/e["v0"]
    
    ax = Axis(f[1,1], ylabel=L"ω_n", xlabel=L"ω", title="J=$J, Dr = $Dr, v0 = $(v0)")

    tau =1/Dr
    eigval_ind=1
    eigval_bin_centers = e["FT_v_projs"]["eigval_bin_centers"]

    w = e["FT_v_projs"]["w"]
    X = zeros(size(e["FT_v_projs"]["X2"]))    

    t = e["t"]
    min_t_ind = e["min_t_ind"]

    a = a_min = sqrt(1 +  (eigval_bin_centers[eigval_ind]/(2 * J) + 1/(2 * tau *J ))^2 ) - (eigval_bin_centers[eigval_ind]/(2 * J) + 1/(2 * tau *J))

    for i in 1:size(X)[1]

        for j in 1:size(X)[2]
            X[i,j] = v0^2  * 2/tau * 2 *pi * w[j]^2/( ( (1/tau + J*a)*eigval_bin_centers[i]-w[j]^2)^2  + w[j]^2 * (eigval_bin_centers[i] - J/a + J*a + 1/tau)^2)*  (length(t)-min_t_ind)/(2*pi)
        end
    end

    heatmap!(ax,w, sqrt.(eigval_bin_centers), log10.(transpose(X)), rasterize=true, colorrange=(-5,5))
    cb = Colorbar(f[1, 2], limits = (-5, 5), colormap = :viridis, flipaxis = false, label = "log10(ω^2 |a_n(ω)|^2)")

    xlims!(0,1)
    ylims!(0,2.5)
    # theory_ABP = v0^2  ./ (2 .+ 2 .* eigval_bin_centers .* tau)  /v0^2

    # scatterlines!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["v0"], colorrange = (0, 0.01 ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")

    # lines!(ax,eigval_bin_centers,theory_ABP, colorrange = (0, maximum(Drs) ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")
    #f[1,2]=Legend(f,ax)
    save("temp.pdf",f)
    append_pdf!( joinpath(figure_save_folder,"wn_w_theory.pdf"), "temp.pdf", cleanup=true)

    display(f)
end


end


















begin
using CairoMakie
CairoMakie.activate!()
f = Figure()
Dr = ensemble_files[1]["Dr"]
ax = Axis(f[1,1], ylabel=L"a", xlabel=L"J", title="Dr = $Dr, v0 = $(v0)")
as = []
a_mins = []
a_mids = []
Js = []
for e in ensemble_files

    vrms = e["vrms"]
    a = vrms/e["v0"]
    a_mid = e["vrms_r"]["vrms_r"][1]

    eigval_bin_centers = e["FT_v_projs"]["eigval_bin_centers"]
    eigval_ind=1
    J = e["J"]
    tau = 1/Dr
    a_min = sqrt(1 +  (eigval_bin_centers[eigval_ind]/(2 * J) + 1/(2 * tau *J ))^2 ) - (eigval_bin_centers[eigval_ind]/(2 * J) + 1/(2 * tau *J))
    push!(a_mins, a_min)
    push!(a_mids, a_mid)
    push!(as, a)
    push!(Js, e["J"])

end

scatterlines!(ax,Js, as)

scatterlines!(ax, Js, a_mins, color="orange", label="a_min")
#scatterlines!(ax, Js, a_mids, color="green", label="a_mid")
ylims!(ax, 0,1)
f[1,2]=Legend(f,ax)
display(f)

end


begin
using CairoMakie
CairoMakie.activate!()
f = Figure()
Dr = ensemble_files[1]["Dr"]
ax = Axis(f[1,1], ylabel=L"a", xlabel=L"J", title="Dr = $Dr, v0 = $(v0)", xscale=log10)
for e in ensemble_files


    eigval_bin_centers = e["FT_v_projs"]["eigval_bin_centers"]
    eigval_ind=1
    J = e["J"]

    Jse = ones(length(eigval_bin_centers)) * J
    tau = 1/Dr
    a_minse = @. sqrt(1 +  (eigval_bin_centers/(2 * J) + 1/(2 * tau *J ))^2 ) - (eigval_bin_centers/(2 * J) + 1/(2 * tau *J))

    scatter!(ax, Jse, a_minse, color=eigval_bin_centers)

end

scatterlines!(ax,Js, as, label="a_num")

scatterlines!(ax, Js, a_mins, color="orange", label="a_min", marker='p')

ylims!(ax, 0,1)
f[1,2]=Legend(f,ax)

save(joinpath(figure_save_folder,"a_vs_J.pdf"), f)

display(f)



end



close.(ensemble_files)