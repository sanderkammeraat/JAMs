

include("../analysis/AnalysisPipeline.jl")

#base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/Nlin_20"

base_folder = "/Users/kammeraat/mounting/data2_kammeraat/sa/statistics/hex_disordered/phi_1.3/Nlin_20"

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

    if Je!=0.03

    scatterlines!(ax, e["vrms_r"]["r_bin_centers"], e["vrms_r"]["vrms_r"]/e["v0"], colorrange = (0, 0.15 ) , color=e["v0"], label="v0 = $(e["v0"]),J = $(e["J"])")
    end

end
f[1,2]=Legend(f,ax)
display(f)

end

begin

f = Figure()

ax = Axis(f[1,1], xlabel=L"λ", ylabel=L"v projs /v0^2", yscale=log10)

for e in ensemble_files


    J= e["J"]
    Dr = e["Dr"]
    v0 = e["v0"]

    if Dr==0.1 && v0<=0.01

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

        scatter!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["v0"], colorrange = (0, 0.01 ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")

        #lines!(ax,eigval_bin_centers,theory_ABP, colorrange = (0, 0.1 ) ,  label="v0 = $(e["v0"]),J = $(e["J"])", color=v0)
        scatterlines!(ax,the_eigvals[select],theory_amin, colorrange = (0, 0.1) ,  label="v0 = $(e["v0"]),J = $(e["J"])", color=v0)
    end
 
end
xlims!(ax, 0,1)
#f[1,2]=Legend(f,ax)
display(f)
end




begin
using CairoMakie
CairoMakie.activate!()
for e in ensemble_files

    f = Figure()

    J= e["J"]
    Dr = e["Dr"]
    v0 = e["v0"]

    ax = Axis(f[1,1], ylabel=L"sqrt(λ)", xlabel=L"ω", title="J=$J, Dr = $Dr, v0 = $(v0)")

    tau =1/Dr

    eigval_bin_centers = e["FT_v_projs"]["eigval_bin_centers"]

    w = e["FT_v_projs"]["w"]

    X = e["FT_v_projs"]["X2"]
    X = e["FT_v_projs"]["X2"]

    heatmap!(ax,w, sqrt.(eigval_bin_centers), log10.(transpose(X)))


    xlims!(0,1)
    # theory_ABP = v0^2  ./ (2 .+ 2 .* eigval_bin_centers .* tau)  /v0^2

    # scatterlines!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["v0"], colorrange = (0, 0.01 ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")

    # lines!(ax,eigval_bin_centers,theory_ABP, colorrange = (0, maximum(Drs) ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")
    #f[1,2]=Legend(f,ax)
    display(f)

end

end












close.(ensemble_files)