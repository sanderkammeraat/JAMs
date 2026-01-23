
using PDFmerger
include("../analysis/AnalysisPipeline.jl")


base_folder = "/Users/kammeraat/mounting/data2_kammeraat/sa/statistics/hex_disordered/phi_1.3/Nlin_20"




figure_save_folder = joinpath(base_folder, "figures_23_01")
mkpath(figure_save_folder)



ef  = findfile(joinpath(base_folder,"ensembles"), "ensemble.h5")


ensemble_files = [ load_file(ef[i]) for i in eachindex(ef)]

Drs = [e["Dr"] for e in ensemble_files]


Js = [e["J"] for e in ensemble_files]

v0s = [e["v0"] for e in ensemble_files] #by color





ensemble_files[1]["eigenmodes"]["eigvals"]["seed_1.h5"]


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

GLMakie.activate!()
f = Figure()

ax = Axis(f[1,1], xlabel=L"$\Delta t$", ylabel=L"$\langle \hat{n}(t+\Delta t)\cdot \hat{n}(t) \rangle$", title="Dr = 0.1")
Jmax = 0.1
for e in ensemble_files

    Dr = e["Dr"]
    J = e["J"]



    if Dr==0.1 && J==0.1

        eigvals = e["eigenmodes"]["eigvals"]["seed_1.h5"]

        #a_min= @.sqrt(1 +  ( eigvals[eigval_ind] / (2 * J) + 1/(2 * tau *J ))^2 ) - ( eigvals[eigval_ind] / (2 * J) + 1/(2 * tau *J ))
        tau = 1/Dr
        a_ABP = sqrt(1/e["Nint"]*sum(1 ./(2 .+ 2*tau .* eigvals)))
        a = a_ABP
        tauJ = 1/(J*a+1/tau)

        det = e["auto_p"]["deltat"]

        the_auto = zeros(length(det))
        #the_auto_s = zeros(length(det))

        for i in eachindex(the_auto)
            ti = det[i]
            #the_auto_s[i] = exp(-ti/tauJ)*tauJ/tau
            the_auto[i] = 2/length(eigvals)*  sum(@.  1/(2*a*tau*(tauJ^2*eigvals^2-1)^2) *exp(-ti*(1/tauJ + eigvals))*( -2*exp(ti/tauJ)*J*tauJ^3*eigvals + exp(ti*eigvals) *tauJ * ( a*(tauJ^2*eigvals^2-1)^2+J*(tauJ + ti+tauJ^2*(tauJ-ti)*eigvals^2)  )   ))
        end
    

    scatter!(ax,det , e["auto_p"]["Cavg"], colorrange = (0, Jmax ) , color=e["J"], label="J = $(e["J"])",alpha=0.3)

    #lines!(ax,det , the_auto, colorrange = (0, Jmax ) , color=e["J"])#, label="J = $(e["J"])")
    end

end

xlims!(ax, 0,100)

lines!(ax, ensemble_files[1]["auto_p"]["deltat"],exp.(-ensemble_files[1]["auto_p"]["deltat"].*0.1), color="red", label="ABP_theory")

ylims!(ax, -0.015,0.2)
f[1,2]=Legend(f,ax)
display(f)
save(joinpath(figure_save_folder,"auto_n_zoom.pdf"), f, backend=CairoMakie)

end
using Roots
function f_a(a,p)

        eigvals =p[1]
        J = p[2]
        tau = p[3] 

        Nint= length(eigvals)/2
        the_eigvals =eigvals[1:end]

        A = Complex.( sqrt.( (1/tau + J * a) .* the_eigvals))

        B = Complex.(@. the_eigvals .+  .-J /a  .+ 1/tau  .+ J*a) 
    
        T1 =@. -(B^2 - 2*A^2)/2
        T2 =@. 1/2*sqrt(B^2 * (B^2 - 4 * A^2)) 
        #theory_integral = @. 2* pi* tauJ * (J*tauJ + a * (1+tauJ*eigvals)^2)/(a*tau*(1+tauJ*eigvals)^3)/4/pi

        theory_integral = @. real( @. 1im * pi* ( 1/ (sqrt(T1+T2)-sqrt(T1-T2)) ) ) *4*pi/tau/(2*pi)^2/2

    return 1/Nint*sum(theory_integral) - a^2
end
using Roots
function f_a_2(a,p)

        eigval =p[1][1]
        J = p[2]
        tau = p[3] 
        num_v = p[4]

        A = Complex.( sqrt.( (1/tau + J * a) .* eigval))

        B = Complex.( eigval .+  .-J /a  .+ 1/tau  .+ J*a) 


    
        T1 =-(B^2 - 2*A^2)/2
        T2 =1/2*sqrt(B^2 * (B^2 - 4 * A^2)) 
        #theory_integral = @. 2* pi* tauJ * (J*tauJ + a * (1+tauJ*eigvals)^2)/(a*tau*(1+tauJ*eigvals)^3)/4/pi

        theory_integral =  real( @. 1im * pi* ( 1/ (sqrt(T1+T2)-sqrt(T1-T2)) ) ) *4*pi/tau/(2*pi)^2/2

        if eigval .+  .-J /a  .+ 1/tau  .+ J*a<=0
            theory_integral*=1000
            
        end

    return theory_integral - num_v
end




display(ensemble_files[1]["t"])
using CairoMakie
GLMakie.activate!()
begin
f = Figure()

ax = Axis(f[1,1], xlabel=L"λ", ylabel=L"v projs /v_0^2",title=L"D_r=0.1", yscale=log10)#, xscale=log10, yscale=log10,title=L"First order in ($J \tau_J$), Dr=0.1")
ϵs = zeros(12)

# J= [0.001, 0.01, 0.02, 0.04, 0.08, 0.10, 0.12, 0.16. 0.2, 0.4, 0.8, 1.0]
ϵs = [0,    0,  0,  0.05, 0.1, 0.2,0.15,0.1, 0.05,0.02,0,0]

ind=1
for e in ensemble_files
    

    
    Dr = e["Dr"]
    v0 = e["v0"]
    J= e["J"]
    if Dr==0.1 && J>0
        display(J)
        
            
        tau =1/Dr

        eigval_bin_centers = e["v_projs_time_avg"]["eigval_bin_centers"]

        eigval_ind = 1 

        

        eigvals = ensemble_files[1]["eigenmodes"]["eigvals"]["seed_2.h5"]

        for key in keys(e["eigenmodes"]["eigvals"])
           # display(e["eigenmodes"]["eigvals"][key][1])
        end

        eigvals = [0.004]
        a_min= @.sqrt(1 +  ( eigvals[eigval_ind] / (2 * J) + 1/(2 * tau *J ))^2 ) - ( eigvals[eigval_ind] / (2 * J) + 1/(2 * tau *J ))
        eigvals =ensemble_files[1]["eigenmodes"]["eigvals"]["seed_2.h5"]# e["eigenmodes"]["eigvals"]["seed_2.h5"]
        a_ABP = sqrt(1/e["Nint"]*sum(1 ./(2 .+ 2*tau .* eigvals))) 

        p= [eigval_bin_centers, J, tau,e["v_projs_time_avg"]["v_projs_time_avg"][1]/e["v0"]^2]
        
        #display(a_num)
        #display(e["vrms"]/v0)
        #a based on loweest mode selection
        #a =copy( e["vrms"]/v0)
        #display(a)
        #a = e["vrms_r"]["vrms_r"][2]/v0
        #a = a_num

        #a = a_alt
        #a=1
        #a = a_num
        #a = a_min

        if a_min<=a_ABP

            a=a_ABP

        else
            a=a_min
        end

        #a+=ϵs[ind]
        ind+=1
        #a_num = find_zero(f_a_2, a, p)
        #display(abs(a_num-e["vrms"]/v0)/a_num*100)
       # a2 = 
        a = a_min
        the_eigvals =eigvals[1:end]
        theory_ABP = v0^2  ./ (2 .+ 2 .* the_eigvals .* tau)  /v0^2
        tauJ = 1/(J*a+1/tau)
        C =@.  1 - 2*(the_eigvals+1/tauJ)*a/J
        A = Complex.( sqrt.( (1/tau + J * a) .* the_eigvals))

        B = Complex.(@. the_eigvals .+  .-J /a  .+ 1/tau  .+ J*a) 
    
        T1 =@. -(B^2 - 2*A^2)/2
        T2 =@. 1/2*sqrt(B^2 * (B^2 - 4 * A^2)) 
        #B=@.B + 0-B[1]

        #integral= @.real( 1im * sqrt(2) *pi/(sqrt(2 *A^2+B*(-B+sqrt(-4 *A^2+B^2)))+sqrt(2 *A^2-B*(B+sqrt(-4* A^2+B^2)))))
        #theory_integral= @.real( 1im  *pi* ( T1+T2 + sqrt(T1-T2)*sqrt(T1+T2))/ (2 *sqrt(T1-T2)*T2) ) ./tau

        theory_integral = @. real( @. 1im * pi* ( 1/ (sqrt(T1+T2)-sqrt(T1-T2)) ) ) *4*pi/tau/(2*pi)^2/2

        
        #theory_integral = @. 2* pi* tauJ * (J*tauJ + a * (1+tauJ*the_eigvals)^2)/(a*tau*(1+tauJ*the_eigvals)^3)

        #theory_integral = @. 4* pi^2 * tauJ/tau * ( 1/(1+tauJ*the_eigvals) + 1* J*tauJ *1/(a*(1+tauJ*the_eigvals)^3))/(2 *pi)^2/2

        #theory_integral_2 =  @. 2* pi* tauJ^2 * (3*J^2 *tauJ +2*J*(1+tauJ *the_eigvals)^3 + 2*(1+tauJ *the_eigvals)^4/tauJ - tauJ*(J+J*tauJ *the_eigvals)^2 /a )/(2*(1+tauJ *the_eigvals)^5*tau)
        
        
        
        #theory_amin = @. pi * ( B - sqrt(B^2 - 4 * A^2))/B/sqrt(-4*A^2 + 2 * B * (B - sqrt(B^2-4*A^2)) ) * 2/tau /4/pi 
        #theory_integral= theory_integral
       # theory_integral_2 = theory_integral_2/4/pi
        #scatter!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2, color=e["J"], colorrange = (0, 1) ,  label="J = $(e["J"])", alpha=0.1)

        #lines!(ax,eigval_bin_centers,theory_ABP, color=e["J"], colorrange = (0, 1) ,  label="J = $(e["J"]) ABP theory ", alpha=0.2)


        #lines!(ax,the_eigvals[select],theory_amin, color=e["J"], colorrange = (0, 1) ,  label="J = $(e["J"]) theory a=a_ABP", linestyle=:dash)

        scatter!(ax,eigval_bin_centers,e["v_projs_time_avg"]["v_projs_time_avg"]/e["v0"]^2,color=log10(e["J"]), colorrange = (-3, 0))# ,  label="J = $(e["J"])",alpha=0.3)

        lines!(ax,the_eigvals,theory_ABP, color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"]) ABP theory ", alpha=0.2)
       
        #scatterlines!(ax,the_eigvals,real.(B), color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"])", linestyle=:dash)

        
        #lines!(ax,the_eigvals,theory_integral, color=log10(e["J"]), colorrange = (-3, 0) , linestyle=:solid, label="J = $(e["J"]), a = $(e["vrms"]/v0)")#,  label="J = $(e["J"])")
        #lines!(ax,the_eigvals,theory_amin, color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"]) theory a=a_ABP", linestyle=:dash)
    end

end

xlims!(ax, 0.002,1)
f[1,2]=Legend(f,ax)
#save(joinpath(figure_save_folder,"v_projs_tau_J_theory_aABP.pdf"), f,backend=CairoMakie)
display(f)#
end


begin #heatmap wn vs w
using CairoMakie
CairoMakie.activate!()
for e in ensemble_files


    f = Figure()

    J= e["J"]
    Dr = e["Dr"]

    if Dr ==0.01 #&& J==0.08
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
    append_pdf!( joinpath(figure_save_folder,"Dre_0_01_wn_w.pdf"), "temp.pdf", cleanup=true)


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
GLMakie.activate!()
CairoMakie.activate!()
begin

f = Figure()

ax = Axis(f[1,1], title="Dr = 0.1", xlabel=L"$ω$", ylabel=L"$\langle |\tilde{p}_x(\omega)|^2\rangle$", xscale=log10, yscale=log10)

for e in ensemble_files

    Dr = e["Dr"]
    J = e["J"]

    Jmax=0.1
    if Dr==0.1 && J==0.1

    eigvals = e["eigenmodes"]["eigvals"]["seed_1.h5"]
    tau=1/Dr
    a_ABP = sqrt(1/e["Nint"]*sum(1 ./(2 .+ 2*tau .* eigvals))) 
    w = e["FT_px"]["w"]
    t = e["t"]
    min_t_ind = e["min_t_ind"]

    
    the_px = zeros(length(w))
    the_px_full = zeros(length(w))
    #a = e["vrms"]/e["v0"]
    a = a_ABP
    tauJ = 1/(a*J+1/tau)
    for i in eachindex(the_px)
        the_px[i] =  1/(length(eigvals)) *sum(  @.  tauJ^2 /tau /(1+tauJ^2 * w[i]^2) * ( 1 +1* tauJ*J * 2*w[i]^2/a /( (1+tauJ^2 * w[i]^2)*(eigvals^2 + w[i]^2) ) ) )* (t[end] - t[min_t_ind])/4
        the_px_full[i] = 1/(length(eigvals)) *sum(  @.         tauJ^2/tau *  (eigvals^2 + w[i]^2)/( (eigvals^2 + w[i]^2)*(1+tauJ^2 * w[i]^2)+w[i]^2*tauJ^2 *J^2/a^2*(1 - 2*(eigvals+1/tauJ)*a/J))         )* (t[end] - t[min_t_ind])/4
    end

    scatter!(ax, e["FT_px"]["w"],e["FT_px"]["X2"], colorrange = (0, Jmax ) , color=(e["J"]), label="J = $(e["J"])",alpha=0.3)

    lines!(ax,w,the_px, colorrange = (0, Jmax ) , color=(e["J"]))#, label="J = $(e["J"])")
    #lines!(ax,w,the_px_full, colorrange = (0, Jmax ) , color=(e["J"]))#, label="J = $(e["J"])")
    end

end
f[1,2]=Legend(f,ax)


xlims!(10^(-3.5), 10^(0.3))

ylims!(10^(4), 10^(4.4))

display(f)
save(joinpath(figure_save_folder,"FT_px_zoom.pdf"), f,backend=CairoMakie)#
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

