
include(joinpath("..","src","Engine.jl"))
include("AnalysisPipeline.jl")
include("AnalysisFunctions.jl")


#base_folder = joinpath(homedir(),"sa","survey","hex_disordered","phi_1","Nlin_4","vary_J_Dr")

#base_folder = joinpath("/Volumes","T7_Shield","sa","single","Dr_0.1","J_0.5_v0_0.3_k_0.4")

#base_folder = joinpath("/Volumes","T7_Shield","sa","single","Dr_0.0001","J_0.5_v0_0.01_k_0.4")
#base_folder=joinpath(homedir(),"test_hdf5", "run")
base_folder = joinpath("/Volumes","T7_Shield","sa","single","Dr_0.01","not_all_1")
#base_folder = "/Volumes/T7_Shield/test_storage/store_vhdf5_v5"

figure_save_folder = mkpath(joinpath(base_folder, "figure_save_folder_23_10"))

#base_folder = joinpath(homedir(),"sa","survey","hex_disordered","phi_1","Nlin_4","vary_J_Dr")
raw_data_base_folder = joinpath(base_folder, "simdata")

raw_data_file_path = joinpath(raw_data_base_folder,"raw_data.jld2")

raw_data_file = jldopen(raw_data_file_path,"r")

JAMS_file =  jldopen(joinpath(raw_data_base_folder,"JAMs_container.jld2"))

frames = raw_data_file["frames"]
system = raw_data_file["system"]
integration_info = raw_data_file["integration_info"]
#frames_support = jldopen(joinpath(raw_data_base_folder,"J_0.1","Dr_0.01","seed_23","ra_raw_data.jld2"),"r")["frames"]

t =  integration_info["save_tax"]
t[2]-t[1]
v0 = frames["1"]["v0"][1]
Dr = frames["1"]["Dr"][1]
J = system["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
k = system["forces"]["external_forces"]["external_harmonic_force"]["k"]
R = frames["1"]["R"]
type = frames["1"]["type"]
Nt = length(t)
#Ignore boundary particles
Np = length(extract_frame_data_for_type("id",1,frames["1"]))


x = zeros(Np, Nt)
y = zeros(Np, Nt)

vx = zeros(Np, Nt)
vy = zeros(Np, Nt)

px = zeros(Np, Nt)
py = zeros(Np, Nt)

qx = zeros(Np, Nt)
qy = zeros(Np, Nt)

@views for i in 1:length(t)
    x[:,i] = extract_frame_data_for_type("x", 1, frames[string(i)])
    y[:,i] = extract_frame_data_for_type("y", 1, frames[string(i)])



    vx[:,i] = extract_frame_data_for_type("vx", 1, frames[string(i)])
    vy[:,i] = extract_frame_data_for_type("vy", 1, frames[string(i)])

    px[:,i] = extract_frame_data_for_type("px", 1, frames[string(i)])
    py[:,i] = extract_frame_data_for_type("py", 1, frames[string(i)])

    qx[:,i] = extract_frame_data_for_type("qx", 1, frames[string(i)])
    qy[:,i] = extract_frame_data_for_type("qy", 1, frames[string(i)])

end


rs =  v0/sqrt(J*k) * sqrt( ( sqrt( 1+(k/(2*J))^2 ) - k/(2*J)) )

r = sqrt.(x.^2 .+ y.^2)


sqrt(J*k*( sqrt(1+(k/(2*J))^2 ) - k/(2*J)))

vs = v0*( ( sqrt( 1+(k/(2*J))^2 ) - k/(2*J)) )
ω_s = sqrt(J*k*vs/v0)


using GLMakie
GLMakie.activate!()

its = 50

# Determine the phase by comparing a unit circle to the unit scale trajectory
phase = angle( x[1, its]/r[1,its] + im * y[1, its]/r[1,its]) - ω_s * t[its]


xs = rs * cos.(ω_s * t .+ phase )

ys = rs * sin.(ω_s * t .+ phase )
#Plot all the particle locations
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"x",ylabel=L"y", aspect=1);


#lines!(ax, x[1,1:1000], y[1,1:1000], color=t[1:1000])

#lines!(ax, [0,0], [0,rs], color="black")

scatter!(ax,x[1,its], y[1,its],label="steady state point")

#scatter!(ax,xs, ys,label="circular orbit positions", color="purple")
lines!(ax, [0,xs[its]],[0, ys[its]], label="circular orbit position at steady state point", color="purple")
later = 2
lines!(ax, [0,xs[its+later]],[0, ys[its+later]], label="circular orbit position $later steps later", color="indigo")
scatter!(ax,x[1,its+later], y[1,its+later],label="steady state point $later steps later", color="pink")

Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)
axislegend()

#ylims!(ax, low=1e-7,high= 1e2)
#save(joinpath(figure_save_folder,"trajectory.pdf"),f)
display(f)
end


begin
f = Figure()
ax = Axis(f[1,1]);

scatter!(ax, t, xs)
scatter!(ax, t, ys)

display(f)
end


begin
f = Figure()
ax = Axis(f[1,1]);

scatter!(ax, t, angle.(  x[1,:] .+ im * y[1,:]))
scatter!(ax, t, angle.(  xs .+ im * ys))
# scatter!(ax , t[its],angle.(  xs .+ im * ys)[its], color="black", markersize=10)
display(f)
end

begin
f = Figure()
ax = Axis(f[1,1]);

scatter!(ax,xs, ys,label="circular orbit positions", color=t)


Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)
axislegend()

#ylims!(ax, low=1e-7,high= 1e2)
#save(joinpath(figure_save_folder,"trajectory.pdf"),f)
display(f)
end


GLMakie.activate!()
function animation(frames,xs, ys, its,rs)

    frame_numbers = its:length(xs)


    t = Observable(0.)

    x = Observable(frames["1"]["x"])
    y = Observable(frames["1"]["y"]) 
    xsO = Observable(xs[1])
    ysO = Observable(ys[1])


    vx = Observable(frames["1"]["vx"])
    vy = Observable(frames["1"]["vy"]) 

    px = Observable(frames["1"]["px"])
    py = Observable(frames["1"]["py"]) 

    R = Observable(frames["1"]["R"]) 

    #Setup figure
    f = Figure(size=(1000,1000));
    ax = Axis(f[1,1], aspect=1,title = @lift("t = $(round($t, digits = 1))"));
    xlims!(ax, (-1,1))

    ylims!(ax, (-1,1))

    c = @lift( angle.($px+1im*$py) )
    scatter!(ax,x,y, markersize=40)

    scatter!(ax,0,0, markersize =2*rs,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))

    scatter!(ax, xsO, ysO, markersize=40, color="red")

    display(f)

    record(f, joinpath(base_folder,"self_alignment_single.mp4"), frame_numbers; framerate=240, visible=true) do i 

        stri = string(i)
        t[] = frames[stri]["t"]

        x[] = frames[stri]["x"]
        y[] = frames[stri]["y"]

        vx[] = frames[stri]["vx"]
        vy[] = frames[stri]["vy"]

        px[] = frames[stri]["px"]
        py[] = frames[stri]["py"]
        xsO[] = xs[i]
        ysO[] = ys[i]

    end
end

animation(frames,xs, ys,its,rs)


α = xs/rs .* x[1,:] + ys/rs .* y[1,:] .- rs


β = -ys/rs .* x[1,:] + xs/rs .* y[1,:] 


dr = r .- rs

dx = x[1,:] .- xs
dy = y[1,:] .- ys
GLMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t",ylabel=L"α", aspect=1);


lines!(ax, t[its:end], α[its:end], color="red", label=L"$α$")

lines!(ax, t[its:end], dr[its:end], color="blue",label=L"$r - r_s $")

lines!(ax, t[its:end], β[its:end], color="green", label=L"$β$")

#lines!(ax, t[its:end], dx[its:end], color="purple")

hlines!(ax, rs,color="black", label= L"$r_s $")

hlines!(ax, -rs,color="black", label= L"$-r_s $")
#lines!(ax, t[its:end], dy[its:end], color="magenta")
#lines!(ax, t[its:end], α[1,its:end], color=t[its:end])


Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)


#ylims!(ax, low=1e-7,high= 1e2)
axislegend()
display(f)
end






GLMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t",ylabel=L"p_x", aspect=1);


lines!(ax, t, px[1,:], color=t)

scatter!(ax, t[its], px[1,its], color="black")



Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)


#ylims!(ax, low=1e-7,high= 1e2)

display(f)
end

#it_max = round(Int64,1500/(t[2] - t[1]))
its=3000
FT_vx = temporal_Fourier_transform(t[2]-t[1],vx,min_t_ind=its)

FT_vy = temporal_Fourier_transform(t[2]-t[1],vy,min_t_ind=its)

using CairoMakie
CairoMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1],yscale=log10,xscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{v_x(t)\}(\omega)|^2");
scatter!(ax, FT_vx["ω"], FT_vx["pavg_X2"], label="vx")

scatter!(ax, FT_vy["ω"], FT_vy["pavg_X2"], color="orange", label="vy")
scatter!(ax,FT_vx["ω_max"] , FT_vx["max_X2"],label="max")

vm = mean(sqrt.(vx.^2 .+ vy.^2)[its:end])  # v0*(sqrt(1+(k/(2*J))^2 ) - k/(2*J))
tau = 1/Dr
a = vm/v0

ω = FT_vx["ω"]

theory_433 =   ω.^2 .* v0^2 * 2/tau * 2 * pi ./ ( ( (1/tau + J * a)*k .- ω.^2).^2 .+ ω.^2 .* (k - J/a + 1/tau + J * a)^2)


theory_433_no_J =   ω.^2 .* v0^2 * 2/tau * 2 * pi ./ ( ( (1/tau + 0*J * a)*k .- ω.^2).^2 .+ ω.^2 .* (k - 0*J/a + 1/tau + 0*J * a)^2)
#lines!(ω, theory_433*(length(t)-its) ,label="theory", linestyle=:dash, color="black")
lines!(ω, theory_433*((length(t)-its)/2/pi),label="theory", linestyle=:solid, color="black")
lines!(ω, theory_433_no_J*((length(t)-its)/2/pi),label="theory, J=0", linestyle=:dash, color="purple")
Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)

#vlines!(ax, sqrt(J*k*( sqrt(1+(k/(2*J))^2 ) - k/(2*J))), label="theory", color="green")
#vlines!(ax, ω_s, label="theory", color="orange")

xlims!(ax, low=1e-2, high=ω[end])

f[1,2]=Legend(f,ax)
#ylims!(ax, low=1e-7,high= 1e2)
#save(joinpath(figure_save_folder,"approximation_analogy_spectrum.pdf"),f)
display(f)
end



#CairoMakie.activate!()
    GLMakie.activate!
begin
f = Figure()
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{x(t)\}(\omega)|^2");
scatter!(ax, FT_x["ω"], FT_x["pavg_X2"])
scatter!(ax,FT_x["ω_max"] , FT_x["max_X2"],label="max")



Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)

vlines!(ax, sqrt(J*k*( sqrt(1+(k/(2*J))^2 ) - k/(2*J))), label="theory", color="green")
vlines!(ax, ω_s, label="theory", color="orange")

vm = mean(sqrt.(vx.^2 .+ vy.^2)[1:19000])  # v0*(sqrt(1+(k/(2*J))^2 ) - k/(2*J))
vlines!(ax, J * sqrt(1 - vm^2/v0^2), label="theory", color="purple")
f[1,2]=Legend(f,ax)
#ylims!(ax, low=1e-7,high= 1e2)
#save(joinpath(figure_save_folder,"single_particle_small_noise_unit_alignment.pdf"),f)
display(f)
end



begin

f = Figure()
ax = Axis(f[1,1]);


Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)



vm =  v0*(sqrt(1+(k/(2*J))^2 ) - k/(2*J))

scatter!(t[1:19000],  sqrt.(vx.^2 .+ vy.^2)[1:19000])

hlines!(ax, vm, label="theory", color="green")
f[1,2]=Legend(f,ax)
#ylims!(ax, low=1e-7,high= 1e2)
#save(joinpath(figure_save_folder,"single_particle_small_noise_unit_alignment.pdf"),f)
display(f)
    
end




FT = temporal_Fourier_transform(t[2]-t[1],px[:,1:it_max],min_t_ind=its)




#CairoMakie.activate!()
GLMakie.activate!
begin
f = Figure()
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{p_x(t)\}(\omega)|^2");
scatter!(ax, FT["ω"], FT["pavg_X2"])
scatter!(ax,FT["ω_max"] , FT["max_X2"],label="max")



Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)


vlines!(ax, sqrt(J*k*( sqrt(1+(k/(2*J))^2 ) - k/(2*J))), label="theory", color="green")
f[1,2]=Legend(f,ax)
#ylims!(ax, low=1e-7,high= 1e2)
#save(joinpath(figure_save_folder,"single_particle_small_noise_unit_alignment.pdf"),f)
display(f)
end
begin
f = Figure()
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{p_x(t)\}(\omega)|^2");
scatter!(ax, FT["ω"], FT["pavg_X2"])
scatter!(ax,FT["ω_max"] , FT["max_X2"],label="max")



Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)


vlines!(ax, sqrt(J*k*( sqrt(1+(k/(2*J))^2 ) - k/(2*J))), label="theory", color="green")
f[1,2]=Legend(f,ax)
ylims!(ax, low=1,high= 1e8)
xlims!(ax, low=0,high= 2)
#save(joinpath(figure_save_folder,"single_particle_small_noise_unit_alignment_zoomed.pdf"),f)
display(f)
end


begin
f = Figure()
ax = Axis(f[1,1], xlabel="t",ylabel=L"< r(t) >");
scatter!(ax,t, r[1,:])
hlines!(ax, rs, color="green")
display(f)
end




begin
f = Figure()
ax = Axis(f[1,1], xlabel="t",ylabel=L"< r(t) >");
scatter!(ax,t, dr[1,:])
hlines!(ax, 0, color="black")
display(f)
#save(joinpath(figure_save_folder,"avg_r_t.pdf"),f)
end

FT_dr = temporal_Fourier_transform(t[2]-t[1],reshape(α, (1, length(α)))[:,1:it_max],min_t_ind=its)
FT_dr_approx = temporal_Fourier_transform(t[2]-t[1],r[:,1:it_max] .- rs,min_t_ind=its)
FT_dr_perp = temporal_Fourier_transform(t[2]-t[1],reshape(β, (1, length(α)))[:,1:it_max],min_t_ind=its)

ω_n = sqrt(J*k*vs/v0)


function theory_spectrum(ω)

    p = v0
    Ω =  sqrt(J*k*vs/v0)

    spectrum=  2*pi* 2 *Dr*( p^2 * Ω ^4 *(J^4 * k^4- 2 * J^2 * k^4 * Ω ^2 + k^2 * (2 *  J^2 + k^2+ω ^2)*  Ω ^4 -2 * k^2 * Ω ^6 + Ω^8 ))/( J^6 * k^6 * (k^2 *  ω ^2 + (ω ^2 - 2*  Ω ^2)^2 )+2 *  J^4 * k^4 * Ω ^4 * ( 3 *  k^2 * ω^2 + (ω^2 - 2 * Ω ^2)^2)  + J^2 *  k^2 * Ω ^4 * (k^4 * ω^4 + Ω^4 *(ω^2 - 2 * Ω^2)^2 + k^2 * (ω^6 - 8 *ω^4 * Ω^2 + 21* ω^2 * Ω^4)) )
return spectrum
end


function theory_perp_spectrum(ω)
    p = v0
    Ω =  sqrt(J*k*vs/v0)
    
    return 2*pi* 2 *Dr* ( p^2 * Ω ^6 * (2 * k^4 * Ω ^2 * (J^2 + 2 * ω ^2 )+ k^2 * Ω ^4 * (2 * J^2+k^2+7* ω^2)+k^4 * (J^4-2 *  J^2 *  ω ^2 + k^2 * ω ^2+ ω^4) + 2 * k^2 *  Ω ^6 + Ω ^8 ))/(J^2 * k^2 * ω ^2 * (J^4 *  k^4 * (k^2 * ω ^2 + (ω^2 - 2 * Ω ^2 )^2 ) + 2 *  J^2 *  k^2 *  Ω ^4 * (3 * k^2 * ω ^2 + (ω ^2 - 2 * Ω ^2 )^2)+Ω^4 * (k^4 * ω ^4 + k^2 * (ω ^6 - 8 *  ω ^4 * Ω ^2 + 21* ω ^2 *  Ω ^4 ) + Ω ^4 * (ω ^2 - 2 * Ω ^2 )^2 )))
end

#(p^2 \[Eta]^2 \[CapitalOmega]^4 (J^4 k^4-2 J^2 k^4 \[CapitalOmega]^2+k^2 (2 J^2+k^2+\[Omega]^2) \[CapitalOmega]^4-2 k^2 \[CapitalOmega]^6+\[CapitalOmega]^8))/(J^6 k^6 (k^2 \[Omega]^2+(\[Omega]^2-2 \[CapitalOmega]^2)^2)+2 J^4 k^4 \[CapitalOmega]^4 (3 k^2 \[Omega]^2+(\[Omega]^2-2 \[CapitalOmega]^2)^2)+J^2 k^2 \[CapitalOmega]^4 (k^4 \[Omega]^4+\[CapitalOmega]^4 (\[Omega]^2-2 \[CapitalOmega]^2)^2+k^2 (\[Omega]^6-8 \[Omega]^4 \[CapitalOmega]^2+21 \[Omega]^2 \[CapitalOmega]^4)))

GLMakie.activate!()

begin
f = Figure()
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{\delta r_{\parallel}(t)\}(\omega)|^2");
scatter!(ax, FT_dr_approx["ω"], FT_dr_approx["pavg_X2"], alpha=0.1,label="single simulation")
#scatter!(ax,FT_dr["ω_max"] , FT_dr["max_X2"],label="max")



Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)

lines!(ax,FT_dr_approx["ω"],  theory_spectrum.(FT_dr_approx["ω"])* (length(t)-its- it_max) , color="green" , label="theory spectrum", linewidth=4)
axislegend()
xlims!(ax, low=1e-5)
#save(joinpath(figure_save_folder,"FT_dr_spectrum_spectrum.pdf"),f)

display(f)
end 
GLMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{\delta r_{\parallel} (t)\}(\omega)|^2");
scatter!(ax, FT_dr_approx["ω"], FT_dr_approx["pavg_X2"], alpha=0.1,label="single simulation")
#scatter!(ax,FT_dr["ω_max"] , FT_dr["max_X2"],label="max")



Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)

lines!(ax,FT_dr_approx["ω"],  theory_spectrum.(FT_dr_approx["ω"])*(length(t)-its- it_max) , color="green" , label="theory spectrum", linewidth=4)
axislegend()
xlims!(ax, low=1e-5, high=4)
#save(joinpath(figure_save_folder,"FT_dr_spectrum_zoomed.pdf"),f)

display(f)
end 




begin
f = Figure()
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{\delta r_{\parallel}(t)\}(\omega)|^2");
scatter!(ax, FT_dr["ω"], FT_dr["pavg_X2"], alpha=0.1,label="single simulation")
#scatter!(ax,FT_dr["ω_max"] , FT_dr["max_X2"],label="max")



Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)

lines!(ax,FT_dr["ω"],  theory_spectrum.(FT_dr["ω"])* (length(t)-its- it_max) , color="green" , label="theory spectrum", linewidth=4)
axislegend()
xlims!(ax, low=1e-5)
#save(joinpath(figure_save_folder,"FT_dr_spectrum_spectrum.pdf"),f)

display(f)
end 
GLMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{\delta r_{\parallel} (t)\}(\omega)|^2");
scatter!(ax, FT_dr["ω"], FT_dr["pavg_X2"], alpha=0.1,label="single simulation")
#scatter!(ax,FT_dr["ω_max"] , FT_dr["max_X2"],label="max")



Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)

lines!(ax,FT_dr["ω"],  theory_spectrum.(FT_dr["ω"])* (length(t)-its- it_max) , color="green" , label="theory spectrum", linewidth=4)
axislegend()
xlims!(ax, low=1e-5, high=4)
#save(joinpath(figure_save_folder,"FT_dr_spectrum_zoomed.pdf"),f)

display(f)
end 

GLMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{\delta r_{\perp}(t)\}(\omega)|^2");
scatter!(ax, FT_dr_perp["ω"], FT_dr_perp["pavg_X2"], alpha=0.1,label="single simulation")
#scatter!(ax,FT_dr["ω_max"] , FT_dr["max_X2"],label="max")



Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)

lines!(ax,FT_dr_perp["ω"],  theory_perp_spectrum.(FT_dr_perp["ω"])* (length(t)-its- it_max) , color="green" , label="theory spectrum", linewidth=4)
axislegend()
xlims!(ax, low=1e-5)
#save(joinpath(figure_save_folder,"FT_dr_spectrum_spectrum.pdf"),f)

display(f)
end 
CairoMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{\delta r_{\perp} (t)\}(\omega)|^2");
scatter!(ax, FT_dr_perp["ω"], FT_dr_perp["pavg_X2"], alpha=0.1,label="single simulation")
#scatter!(ax,FT_dr["ω_max"] , FT_dr["max_X2"],label="max")



Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)

lines!(ax,FT_dr_perp["ω"],  theory_perp_spectrum.(FT_dr_perp["ω"])* (length(t)-its- it_max) , color="green" , label="theory spectrum", linewidth=4)
axislegend()
xlims!(ax, low=1e-5, high=4)
#save(joinpath(figure_save_folder,"FT_dr_spectrum_zoomed.pdf"),f)

display(f)
end 




xg = collect(range(0,10,1000))

yg = ones((1,length(xg)))
w= 30
yg[1,:]=cos.(w.* xg) 

FT_yg = temporal_Fourier_transform(xg[2]-xg[1], yg)


begin
f = Figure()
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{p_x(t)\}(\omega)|^2");

scatter!(ax, FT_yg["ω"], FT_yg["pavg_X2"])

hlines!(ax, (yg[1]*2*pi*1/2)^2 *(length(xg)/2/pi)^2   , color="red")
vlines!(ax,w)
display(f)

end 