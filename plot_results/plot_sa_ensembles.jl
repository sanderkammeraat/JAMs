

include("../analysis/AnalysisPipeline.jl")

base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/Nlin_20"

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

    if Dr==0.01

        tau =1/Dr

        eigval_bin_centers = e["v_projs_time_avg"]["eigval_bin_centers"]

        theory_ABP = v0^2  ./ (2 .+ 2 .* eigval_bin_centers .* tau)  /v0^2

        scatterlines!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["v0"], colorrange = (0, 0.01 ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")

        lines!(ax,eigval_bin_centers,theory_ABP, colorrange = (0, maximum(Drs) ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")

    end


end
f[1,2]=Legend(f,ax)
display(f)

end











close.(ensemble_files)