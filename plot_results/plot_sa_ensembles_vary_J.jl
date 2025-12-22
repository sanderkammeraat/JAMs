
using PDFmerger
include("../analysis/AnalysisPipeline.jl")

#base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/Nlin_20"

base_folder = "/Users/kammeraat/mounting/data2_kammeraat/sa/statistics/hex_disordered/phi_1.3/Nlin_20"




figure_save_folder = joinpath(base_folder, "figures_19_12")
mkpath(figure_save_folder)



ef  = findfile(base_folder, "ensemble.h5")


ensemble_files = [ load_file(ef[i]) for i in eachindex(ef)]

Drs = [e["Dr"] for e in ensemble_files]


Js = [e["J"] for e in ensemble_files]

v0s = [e["v0"] for e in ensemble_files] #by color








GLMakie.activate!()
begin

f = Figure()

ax = Axis(f[1,1], xlabel=L"r", ylabel=L"vrms(r)/v0")

for e in ensemble_files

    Dre = e["Dr"]

    if Dre==0.1
    Je = e["J"]

    scatterlines!(ax, e["vrms_r"]["r_bin_centers"], e["vrms_r"]["vrms_r"]/e["v0"], colorrange = (0, 1 ) , color=e["J"], label="J = $(e["J"])")
    end

end
ylims!(ax,0,1.5)
f[1,2]=Legend(f,ax)
display(f)

end


begin
using CairoMakie
CairoMakie.activate!()


f = Figure()

ax = Axis(f[1,1], xlabel=L"$\Delta t$", ylabel=L"$\langle \hat{n}(t+\Delta t)\cdot \hat{n}(t) \rangle$", title="Dr = 0.1")
Jmax = 0.1
for e in ensemble_files

    Dre = e["Dr"]
    Je = e["J"]



    if Dre==0.1 && Je<=Jmax
    

    scatterlines!(ax, e["auto_p"]["deltat"], e["auto_p"]["Cavg"], colorrange = (0, Jmax ) , color=e["J"], label="J = $(e["J"])")
    end

end

xlims!(ax, 0,200)

lines!(ax, ensemble_files[1]["auto_p"]["deltat"],exp.(-ensemble_files[1]["auto_p"]["deltat"].*0.1), color="red", label="ABP_theory")

ylims!(ax, -0.02,0.02)
f[1,2]=Legend(f,ax)
display(f)
save(joinpath(figure_save_folder,"auto_n_zoom.pdf"), f)

end


begin


f = Figure()

ax = Axis(f[1,1], xlabel=L"λ", ylabel=L"v projs /v_0^2")

for e in ensemble_files


    
    Dr = e["Dr"]
    v0 = e["v0"]

    if Dr==0.1# && J==0.0
        J= e["J"]
            
        tau =1/Dr

        eigval_bin_centers = e["v_projs_time_avg"]["eigval_bin_centers"]

        eigval_ind = 1 

        theory_ABP = v0^2  ./ (2 .+ 2 .* eigval_bin_centers .* tau)  /v0^2

        eigvals = e["eigenmodes"]["eigvals"]["seed_1.h5"]


        a= @.sqrt(1 +  ( eigvals[eigval_ind] / (2 * J) + 1/(2 * tau *J ))^2 ) - ( eigvals[eigval_ind] / (2 * J) + 1/(2 * tau *J ))
        a_ABP = sqrt(1/e["Nint"]*sum(1 ./(2 .+ 2*tau .* eigvals))) 
        #a based on loweest mode selection
        #a = e["vrms"]/v0
        #display(a)
        #a = e["vrms"]/v0
        a = a_ABP

        the_eigvals = eigvals[1:end]

        A =@. Complex( sqrt.( (1/tau + J * a) .* the_eigvals))

        B =@. Complex(the_eigvals .+  .-J /a  .+ 1/tau  .+ J*a)
        #B=@.B + 0-B[1]

        I= @.real( 1im * sqrt(2) *pi/(sqrt(2 *A^2+B*(-B+sqrt(-4 *A^2+B^2)))+sqrt(2 *A^2-B*(B+sqrt(-4* A^2+B^2)))))
        #T1 = -(B^2 - 2*A^2)/2
        #T2 = 1/2*sqrt(Complex(B^2 * (B^2 - 4 * A^2))) 

        
        theory_amin = @. pi * ( B - sqrt(B^2 - 4 * A^2))/B/sqrt(-4*A^2 + 2 * B * (B - sqrt(B^2-4*A^2)) ) * 2/tau /4/pi 
        theory_I= I* 2/tau /4/pi 
        #scatter!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["J"], colorrange = (0, 1) ,  label="J = $(e["J"])", alpha=0.1)

        #lines!(ax,eigval_bin_centers,theory_ABP, color=e["J"], colorrange = (0, 1) ,  label="J = $(e["J"]) ABP theory ", alpha=0.2)


        #lines!(ax,the_eigvals[select],theory_amin, color=e["J"], colorrange = (0, 1) ,  label="J = $(e["J"]) theory a=a_ABP", linestyle=:dash)

        #scatterlines!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"])")

        #lines!(ax,eigval_bin_centers,theory_ABP, color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"]) ABP theory ", alpha=0.2)
        scatterlines!(ax,the_eigvals,real.(B), color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"])", linestyle=:dash)
        #lines!(ax,the_eigvals,theory_I, color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"])", linestyle=:dash)
        #lines!(ax,the_eigvals,theory_amin, color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"]) theory a=a_ABP", linestyle=:dash)
    end

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

    if Dr ==0.1 #&& J==0.08
    v0 = e["v0"]

    ax = Axis(f[1,1], ylabel=L"ω_n", xlabel=L"ω", title="J=$J, Dr = $Dr, v0 = $(v0)")

    tau =1/Dr

    eigval_bin_centers = e["FT_v_projs"]["eigval_bin_centers"]

    w = e["FT_v_projs"]["w"]

    X = e["FT_v_projs"]["X2"]

    heatmap!(ax,w, sqrt.(eigval_bin_centers), log10.(transpose(X)), rasterize=true, colorrange=(-2,4), colormap=:gist_rainbow)
    cb = Colorbar(f[1, 2], limits = (-2, 4), colormap = :gist_rainbow, flipaxis = false, label = "log10(ω^2 |a_n(ω)|^2)")


    xlims!(0,0.5)
    ylims!(0,1)
    # theory_ABP = v0^2  ./ (2 .+ 2 .* eigval_bin_centers .* tau)  /v0^2

    # scatterlines!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["v0"], colorrange = (0, 0.01 ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")

    # lines!(ax,eigval_bin_centers,theory_ABP, colorrange = (0, maximum(Drs) ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")
    #f[1,2]=Legend(f,ax)
    save("temp.pdf",f)
    append_pdf!( joinpath(figure_save_folder,"Dre_0_1_wn_w.pdf"), "temp.pdf", cleanup=true)


    display(f)
    end
end
end

e1 = ensemble_files[1]

eigvals = e1["eigenmodes"]["seed_1.h5"]
tau = 1/e1["Dr"]
a =  sqrt(1/e1["Nint"]*sum(1 ./(2 .+ 2*tau .* eigvals)) ) #vrms/e["v0"]

begin #heatmap wn vs w THEORY
using CairoMakie
CairoMakie.activate!()
for e in ensemble_files

    f = Figure()

    J= e["J"]
    Dr = e["Dr"]
    v0 = e["v0"]
    tau =1/Dr
    vrms = e["vrms"]
    
    if Dr ==0.1# && J==0.0

    
    ax = Axis(f[1,1], ylabel=L"ω_n", xlabel=L"ω", title="J=$J, Dr = $Dr, v0 = $(v0)")

    
    eigval_ind=1
    eigval_bin_centers = e["FT_v_projs"]["eigval_bin_centers"]

    eigvals = e["eigenmodes"]["eigvals"]["seed_1.h5"]
    


    w = e["FT_v_projs"]["w"]
    X = zeros(size(e["FT_v_projs"]["X2"]))    

    t = e["t"]
    min_t_ind = e["min_t_ind"]

    #a = a_min = sqrt(1 +  (eigval_bin_centers[eigval_ind]/(2 * J) + 1/(2 * tau *J ))^2 ) - (eigval_bin_centers[eigval_ind]/(2 * J) + 1/(2 * tau *J))
    a = a_ABP = sqrt(1/e["Nint"]*sum(1 ./(2 .+ 2*tau .* eigvals))) #
    #a = vrms/v0

    for i in 1:size(X)[1]

        for j in 1:size(X)[2]
            #X[i,j] = v0^2  * 2/tau * 2 *pi * w[j]^2/( ( (1/tau + J*a)*eigval_bin_centers[i]-w[j]^2)^2  + w[j]^2 * (eigval_bin_centers[i] - J/a + J*a + 1/tau)^2) * (t[end] - t[min_t_ind])/(2*pi)/(2*pi)
            X[i,j] = v0^2  * 2*tau * 2 *pi * w[j]^2*(1)/( (eigvals[i]^2 + w[j]^2)*(1+tau^2 * w[j]^2)-2*J*tau/a*((1+eigvals[i]*tau)*w[j]^2 - eigvals[i]^2*a^2 ) ) * (t[end] - t[min_t_ind])/(2*pi)/(2*pi)

        end
    end

    heatmap!(ax,w, sqrt.(eigval_bin_centers), log10.(transpose(X)), rasterize=true, colorrange=(-2,4),colormap=:gist_rainbow)
    cb = Colorbar(f[1, 2], limits = (-2, 4), colormap = :gist_rainbow, flipaxis = false, label = "log10(ω^2 |a_n(ω)|^2)")

    xlims!(0,0.5)
    ylims!(0,1.)
    # theory_ABP = v0^2  ./ (2 .+ 2 .* eigval_bin_centers .* tau)  /v0^2

    # scatterlines!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["v0"], colorrange = (0, 0.01 ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")

    # lines!(ax,eigval_bin_centers,theory_ABP, colorrange = (0, maximum(Drs) ) ,  label="v0 = $(e["v0"]),J = $(e["J"])")
    #f[1,2]=Legend(f,ax)
    #save("temp.pdf",f)
    #append_pdf!( joinpath(figure_save_folder,"Dre_0.1_wn_w_theory_a_ABP.pdf"), "temp.pdf", cleanup=true)

    display(f)
    end
end


end






















begin

    e = ensemble_files[2]

    f = Figure()

    J= e["J"]
    Dr = e["Dr"]
    v0 = e["v0"]

    ax = Axis(f[1,1], ylabel=L"ω_n", xlabel=L"ω", title="J=$J, Dr = $Dr, v0 = $(v0)")

    tau =1/Dr

    eigval_bin_centers = e["FT_v_projs"]["eigval_bin_centers"]

    w = e["FT_v_projs"]["w"]

    X = e["FT_v_projs"]["X2"]

    heatmap!(ax,w, sqrt.(eigval_bin_centers), log10.(transpose(X)), rasterize=true, colorrange=(-2,5), colormap=:gist_rainbow)
    cb = Colorbar(f[1, 2], limits = (-2, 5), colormap = :gist_rainbow, flipaxis = false, label = "log10(ω^2 |a_n(ω)|^2)")
    
end



begin
using CairoMakie
CairoMakie.activate!()
for Drplot in [0.1,0.01]
    f = Figure()
    v0=ensemble_files[1]["v0"]

    ax = Axis(f[1,1], ylabel=L"a", xlabel=L"J", title="Dr = $Drplot, v0 = $(v0)", xscale=log10)
    as = []
    a_mins = []
    a_mids = []
    eigvals = ensemble_files[2]["eigenmodes"]["eigvals"]["seed_1.h5"]
    tau = 1/Drplot
    a_ABP = sqrt(1/ensemble_files[1]["Nint"]*sum(1 ./(2 .+ 2*tau .* eigvals))) 
    Js = []
    for e in ensemble_files

        Dr = e["Dr"]

        if Dr==Drplot

        vrms = e["vrms"]
        a = vrms/e["v0"]


        a_mid = e["vrms_r"]["vrms_r"][1]

        eigval_bin_centers = e["FT_v_projs"]["eigval_bin_centers"]
        eigval_ind=1
        J = e["J"]
        tau = 1/Dr
        a_min = sqrt(1 +  (eigvals[eigval_ind]/(2 * J) + 1/(2 * tau *J ))^2 ) - (eigvals[eigval_ind]/(2 * J) + 1/(2 * tau *J))
        push!(a_mins, a_min)
        push!(a_mids, a_mid)
        push!(as, a)
        push!(Js, e["J"])

        Jse = ones(length(eigvals)) * J
        tau = 1/Dr
        a_minse = @. sqrt(1 +  (eigvals/(2 * J) + 1/(2 * tau *J ))^2 ) - (eigvals/(2 * J) + 1/(2 * tau *J))

        #a_minse_2 = @. J/(1/tau + eigvals)

        scatter!(ax, Jse, a_minse, color=eigvals)


        end

    end

    scatterlines!(ax,Js, as, label="a_num")

    scatterlines!(ax, Js, a_mins, color="orange", label="a_min", marker='p')
    hlines!(ax, a_ABP, color="red", label="a_ABP")

    annotation!(ax, (0.04,0), text="0.04")
    annotation!(ax, (0.08,0), text="0.08")
    #scatterlines!(ax, Js, a_mids, color="green", label="a_mid")
    ylims!(ax, 0,1)
    f[1,2]=Legend(f,ax)
    display(f)

    save(joinpath(figure_save_folder,"a_vs_J_Dr_$(Drplot).pdf"), f)
end
end


begin
using CairoMakie
CairoMakie.activate!()
f = Figure()
Dr = ensemble_files[1]["Dr"]
v0=0.01
ax = Axis(f[1,1], ylabel=L"a", xlabel=L"J", title="Dr = $Dr, v0 = $(v0)", xscale=log10)
for e in ensemble_files


    eigval_bin_centers = e["FT_v_projs"]["eigval_bin_centers"]
    eigval_ind=1
    J = e["J"]

    if J==0.001 || J==0.0
    Jse = ones(length(eigval_bin_centers)) * J
    tau = 1/Dr
    a_minse = @. sqrt(1 +  (eigval_bin_centers/(2 * J) + 1/(2 * tau *J ))^2 ) - (eigval_bin_centers/(2 * J) + 1/(2 * tau *J))
    a_minse2 = @. 
    scatter!(ax, Jse, a_minse, color=eigval_bin_centers)
    end
end

scatterlines!(ax,Js, as, label="a_num")

scatterlines!(ax, Js, a_mins, color="orange", label="a_min", marker='p')

ylims!(ax, 0,1)
f[1,2]=Legend(f,ax)

save(joinpath(figure_save_folder,"a_vs_J.pdf"), f)

display(f)



end





begin
using GLMakie
GLMakie.activate!()
for e in ensemble_files


    J= e["J"]
    Dr = e["Dr"]
    v0 = e["v0"]
    tau =1/Dr
    vrms = e["vrms"]
    
    if Dr ==0.1 && J==0.0
        f = Figure()
        
        ax = Axis(f[1,1], ylabel=L"ω_n", xlabel=L"ω", title="J=$J, Dr = $Dr, v0 = $(v0)",yscale=log10)

        
        eigval_ind=1
        eigval_bin_centers = e["FT_v_projs"]["eigval_bin_centers"]

        eigvals = e["eigenmodes"]["eigvals"]["seed_1.h5"]
        


        w = e["FT_v_projs"]["w"]
        X = zeros(size(e["FT_v_projs"]["X2"]))    

        t = e["t"]
        min_t_ind = e["min_t_ind"]

        a = e["vrms"]/v0
        display(a)
        i=1

        theory = @. v0^2  * 2/tau * 2 *pi * w^2/( ( (1/tau + J*a)*eigval_bin_centers[i]-w^2)^2  + w^2 * (eigval_bin_centers[i] - J/a + J*a + 1/tau)^2).* (t[end] - t[min_t_ind])/(2*pi)^2#
        #theory = theory .* (t[end] - t[min_t_ind])/(2*pi)^2
        #theory = theory .* (length(t) - min_t_ind)/(2*pi)/2/pi
        display(t[end] - t[min_t_ind])
        display(length(t) - min_t_ind)
        numerics = e["FT_v_projs"]["X2"][i,:]
        
        scatter!(ax, w, numerics)
        lines!(ax, w, theory, color="black")
    end
end
display(f)
end



using CairoMakie
CairoMakie.activate!()
begin

f = Figure()

ax = Axis(f[1,1], title="Dr = 0.1", xlabel=L"$ω$", ylabel=L"$\langle |\tilde{p}_x(\omega)|^2\rangle$", xscale=log10, yscale=log10)

for e in ensemble_files

    Dre = e["Dr"]
    Je = e["J"]

    Jmax=0.1
    if Dre==0.1 && Je<=Jmax
    

    scatterlines!(ax, e["FT_px"]["w"],e["FT_px"]["X2"][1,:], colorrange = (0, Jmax ) , color=(e["J"]), label="J = $(e["J"])")
    end

end
f[1,2]=Legend(f,ax)


xlims!(10^(-3.5), 10^(0.3))

ylims!(10^(4), 10^(4.4))

display(f)
save(joinpath(figure_save_folder,"FT_px_zoom.pdf"), f)

end


using Integrals

function theory_spectrum_n(w,pars)


        eigval_n = pars[1]

        J = pars[2]
        a = pars[3]
        v0 = pars[4]
        tau = pars[5]

    return 2/tau * 2 *pi * w^2/( ( (1/tau + J*a)*eigval_n-w^2)^2  + w^2 * (eigval_n - J/a + J*a + 1/tau)^2)/(2*pi)/(2*pi)
end

begin 
f = Figure()
ax = Axis(f[1,1], yscale=log10)

for e in ensemble_files
    eigvals = e["eigenmodes"]["eigvals"]["seed_1.h5"]
    eigval_ind = 1
    t = e["t"]
    Dr = e["Dr"]
    tau = 1/Dr
    J = e["J"]
    display(J)





    if J>0. && Dr>=0.1
        v0 = e["v0"]
        a = sqrt(1/e["Nint"]*sum(1 ./(2 .+ 2*tau .* eigvals))) # # (sqrt( 1 +  ( eigvals[eigval_ind]/(2 * J) + 1/(2 * tau *J ) )^2 ) - (eigvals[eigval_ind]/(2 * J) + 1/(2 * tau *J)))*1.01
        display(a)
        min_t_ind = e["min_t_ind"]

        sols = []
        

        for eigval in eigvals[2:end]

            domain = (-Inf, Inf)

            pars = [eigval, J, a, v0, tau]
            prob = IntegralProblem(theory_spectrum_n, domain, pars)
            push!(sols,solve(prob, QuadGKJL()).u)
        end

        eigval_bin_centers = e["v_projs_time_avg"]["eigval_bin_centers"]

        scatterlines!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["J"], colorrange = (0, 1) ,  label="J = $(e["J"])")
        lines!(ax, eigvals[2:end], sols, color=e["J"], colorrange = (0, 1))



    end

end
xlims!(ax, 1e-3,2)
f[1,2]=Legend(f,ax)
display(f)

end
close.(ensemble_files)

prob = IntegralProblem(theory_spectrum_n, domain, pars)

