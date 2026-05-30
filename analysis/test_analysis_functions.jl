
begin
include(joinpath("..","src","Engine.jl"))
include("AnalysisPipeline.jl")
include("AnalysisFunctions.jl")
using ProgressMeter
using CairoMakie
CairoMakie.activate!()

#base_folder = joinpath(homedir(),"sa","survey","hex_disordered","phi_1","Nlin_4","vary_J_Dr")

base_folder = joinpath("/Volumes","T7_Shield","sa","survey","hex_disordered", "phi_1", "Nlin_20", "vary_J_Dr")

#base_folder = joinpath("/Volumes","T7_Shield","sa","survey","hex_disordered", "phi_1", "Nlin_50")


#base_folder =  joinpath(homedir(),"mounting", "alicedata1_kammeraatsc1",  "JuliaOnALICE","hex_dis","Nlin_100")

figure_save_folder = mkpath(joinpath(base_folder,"exploratory_figures_14_11"))

#base_folder = joinpath(homedir(),"sa","survey","hex_disordered","phi_1","Nlin_4","vary_J_Dr")
#raw_data_base_folder = joinpath(base_folder, "simdata")

raw_data_base_folder = joinpath(base_folder, "simdata","J_1.0","Dr_0.01", "seed_53")

#Make tree to navigate simulation data folder structure
#tree = construct_folder_tree_param_param_seed(raw_data_base_folder)

#raw_data_file_path = joinpath(raw_data_base_folder,"J_0.1","Dr_0.01","seed_1","sa_raw_data.h5")

raw_data_file_path = joinpath(raw_data_base_folder,"sa_raw_data.jld2")


raw_data_file = jldopen(raw_data_file_path, "r")

frames = raw_data_file["frames"]
system = raw_data_file["system"]

print(frames)

#frames_support = jldopen(joinpath(raw_data_base_folder,"J_0.1","Dr_0.01","seed_1","ra_raw_data.h5"),"r")["frames"]

#For now use last frame of first relaxation step
frames_support = jldopen(joinpath(raw_data_base_folder,"rx_raw_data.jld2"),"r")["frames"]

t =  raw_data_file["integration_info"]["save_tax"]
v0 = frames["1"]["v0"][1]
Dr = frames["1"]["Dr"][1]
J = system["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]
k = system["forces"]["pair_forces"]["soft_disk_force"]["karray"]
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


@views for i in 1:Nt
    x[:,i] .= extract_frame_data_for_type("x", 1, frames[string(i)])
    y[:,i] .= extract_frame_data_for_type("y", 1, frames[string(i)])



    vx[:,i] .= extract_frame_data_for_type("vx", 1, frames[string(i)])
    vy[:,i] .= extract_frame_data_for_type("vy", 1, frames[string(i)])

    px[:,i] .= extract_frame_data_for_type("px", 1, frames[string(i)])
    py[:,i] .= extract_frame_data_for_type("py", 1, frames[string(i)])

    qx[:,i] .= extract_frame_data_for_type("qx", 1, frames[string(i)])
    qy[:,i] .= extract_frame_data_for_type("qy", 1, frames[string(i)])

end


x0 = frames_support[string(length(frames_support))]["x"]
y0 = frames_support[string(length(frames_support))]["y"]

x0int = extract_frame_data_for_type("x",1,frames_support[string(length(frames_support))])
y0int = extract_frame_data_for_type("y",1,frames_support[string(length(frames_support))])

end

include("AnalysisFunctions.jl")
begin
#With boundary
display("Constructing D_wb")
D_wb = construct_D(x0, y0, k, R, type)


#@profview_allocs D_wb = construct_D(x0, y0, k, R, type)
#Boundary particles are stored at the end
interior_indmax = sum(type.==1)

display("Slicing out D")
D = Symmetric(D_wb[1:2*interior_indmax,1:2*interior_indmax])

display("Diagonalizing D")
eigenmodes = diagonalize_D(D)

end


function get_tag()


    k =system["forces"]["pair_forces"]["soft_disk_force"]["karray"]

    #Interior particles
    Nint = sum(type .== 1)
    ϕ = 1
    #Check if all radii are the same, if so do, else , hardcoded 0.15, because I did not store the polydispersity of the initial conditions in the  analysis file
    poly =  all( R .== R[1]) ? 0. : 0.15


    tag = Dict("ϕ"=>ϕ, "v0"=> v0, "Nint"=> Nint, "poly"=>poly, "k"=>k, "Dr"=>Dr, "J"=>J)

    return tag

end
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L"v/v_0")
for i in 1:100:size(vx)[1]
    lines!(ax, t[500:end], sqrt.(vx.^2 + vy.^2)[i,500:end]/v0)
end
#modenumbers = range(1,size(v_projs)[1])
#heatmap!(ax, sqrt.(eigenmodes["eigvals"]),FT_v_projs["ω"], log10.(FT_v_projs["Xf2"]))
tag = get_tag()
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
#f[1,2]=Legend(f,ax)
save(joinpath(figure_save_folder,"mean_v_prediction_constant_speed_particles.pdf"),f)
display(f)#

end

using LinearAlgebra


P = D - diagm(1=>diag(D, 1)) - diagm(-1=>diag(D, -1)) - diagm(diag(D, 0))

K = D - P

GLMakie.activate!()
begin
    heatmap(P)
end

begin
    heatmap(D)
end

begin
    heatmap(K)
end


i=Nt

ri = sqrt.( x[:,i] .^2 + y[:,i].^2 )

xyi = collect(Iterators.flatten(zip(x[:,i] - x0int, y[:,i] - y0int)))

vxyi = collect(Iterators.flatten(zip(vx[:,i], vy[:,i])))

PdR = P * xyi

PdRdot = P * vxyi

normPdR = norm(PdR)
normPdRdot = norm(PdRdot)

normdR = norm(xyi)


vi = sqrt.(vxyi .* vxyi)

C = diagm(vi)

ns = J/v0.* C * PdR + PdRdot


norm(ns)

GLMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1])
scatter!(ax,J/v0 .* C * K*xyi, color="blue")

scatter!(ax,ns, color="red")    
scatter!(ax,J/v0.* C * PdR, color="green")
scatter!(ax,PdRdot, color="purple")
scatter!(ax,K * vxyi, color="pink")

scatter!(ax,K * vxyi - J*v0*inv(C) * vxyi + J/v0 * C * vxyi, color="orange")



hlines!(ax,norm(J/v0 .* C * K*xyi), color="blue")
hlines!(ax,norm(ns), color="red")
hlines!(ax,norm(J/v0.* C * PdR), color="green")
hlines!(ax,norm(PdRdot), color="purple")

display(f)
end

GLMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1])
scatter!(ax,mean(sqrt.(vx.^2 + vy.^2), dims=2)[:,1]/v0, color=ri)
scatter!(ax, (sqrt.( 1 .+ (diag(K)[1:2:end]/(2*J)) .^2 ) - diag(K)[1:2:end]/(2*J)) )



display(f)
end
begin
f = Figure()
ax = Axis(f[1,1])
scatter!(ax,mean(sqrt.((x .- x0int).^2 + (y .- y0int).^2), dims=2)[:,1], color=ri)
scatter!(ax, v0* sqrt.((sqrt.( 1 .+ (diag(K)[1:2:end]/(2*J)) .^2 ) - diag(K)[1:2:end]/(2*J)))/sqrt.(J .* diag(K)[1:2:end] ) )
#


display(f)
end






CairoMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L"v/v_0")

scatter!(ax, t, mean(sqrt.(vx.^2 + vy.^2),dims=1)[1,:]/v0,label="average speed")

scatter!(ax, t, sqrt.( mean((vx.^2 + vy.^2),dims=1)[1,:] )/v0,label="rms speed")
vmean = mean(sqrt.( mean((vx.^2 + vy.^2), dims=1)[1,500:end] ))
interval = 1
# for (n,eigval) in pairs(eigenmodes["eigvals"][1:interval:10])
#     hlines!(ax, ( sqrt( 1 + (eigval/(2*J))^2 ) - eigval/(2*J) ),color="black", label="mode number $(1+(n-1)*interval)")
# end

hlines!(ax, vmean/v0, color="black")
#xlims!(ax, 5000,t[end])
ylims!(0,1)
#modenumbers = range(1,size(v_projs)[1])
#heatmap!(ax, sqrt.(eigenmodes["eigvals"]),FT_v_projs["ω"], log10.(FT_v_projs["Xf2"]))
tag = get_tag()
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

f[1,2]=Legend(f,ax)
save(joinpath(figure_save_folder,"mean_v_prediction_constant_speed.pdf"),f)
display(f)#
end


dt = t[2] - t[1]


FT_speed = temporal_Fourier_transform(dt,  sqrt.( mean((vx.^2 + vy.^2),dims=1)), min_t_ind = 40000, output_not_avg=true)


GLMakie.activate!()
begin
f = Figure()#
ax = Axis(f[1,1], xlabel="ω",ylabel=L"|  \mathcal{F}v(t)\}(\omega)|^2", yscale=log10);
scatter!(ax, FT_speed["ω"][2:end], FT_speed["Xf2"][1,2:end], label="FT{v}")
#scatter!(ax,FT_speed["ω_max"] , FT_speed["max_X2"],label="max")


#vlines!(ax, J*sqrt(1 - vmean^2/v0^2),color="red")
Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)
#xlims!(ax, 0,0.5)

#vlines!(ax, sqrt(J*k*( sqrt(1+(k/(2*J))^2 ) - k/(2*J))), label="theory", color="green")
f[1,2]=Legend(f,ax)
#ylims!(ax, low=1e-7,high= 1e2)
#save("single_particle_small_noise_unit_alignment.pdf",f)
display(f)
end


dx = x .- x0int
dy = y .- y0int
d_projs = project_on_eigvecs(eigenmodes["eigvecs"], dx, dy)

v_projs = project_on_eigvecs(eigenmodes["eigvecs"], vx,vy)



p_projs = project_on_eigvecs(eigenmodes["eigvecs"], px,py)

GLMakie.activate!()
#CairoMakie.activate!()
begin
f = Figure()
ωs = sqrt.(eigenmodes["eigvals"])
ax = Axis(f[1,1], xlabel=L"ω_n", ylabel=L"vproj^2", yscale=log10, xscale=log10)
                        
tau =1/Dr
theory_ABP = v0^2  ./ (2 .+ 2 .* ωs.^2 .* tau)



prefactor = 1/4/pi

eigval_ind=1
a_min = sqrt(1 +  (eigenmodes["eigvals"][eigval_ind]/(2 * J) + 1/(2 * tau *J ))^2 ) - (eigenmodes["eigvals"][eigval_ind]/(2 * J) + 1/(2 * tau *J))

#a based on loweest mode selection
a = a_min


A =  sqrt.((1/tau + J * a) .* eigenmodes["eigvals"])

B = eigenmodes["eigvals"] .+ ( -J/a  + 1/tau  + J*a)

display(B[1])
display(4*A[1]- B[1]^2)
display(B[1]^2)


select = B .> 2* A

A =   A[select]
B = B[select]
theory_amin = @. pi * (B - sqrt(B^2 - 4 * A^2))/B/sqrt(-4*A^2 + 2 * B * (B - sqrt(B^2-4*A^2))) * v0^2 *2/tau*prefactor

numerics =  mean(v_projs[:,1:end].^2, dims=2)[:,1]
scatterlines!(ax,ωs[1:end], numerics[1:end],  label="numerics", linewidth=0.5, alpha=0.1 )

lines!(ax, ωs[select], theory_amin, alpha=1, color="black", label="lowest mode -> a")
lines!(ax, ωs, theory_ABP, linewidth=5, color="red", label="theory ABP")
#a based on empirically most excited mode selection
eigval_ind=findmax(numerics)[2]
a_min = sqrt(1 +  (eigenmodes["eigvals"][eigval_ind]/(2 * J) + 1/(2 * tau *J ))^2 ) - (eigenmodes["eigvals"][eigval_ind]/(2 * J) + 1/(2 * tau *J))
a = a_min
A = sqrt.( (1/tau + J * a) .* eigenmodes["eigvals"])

B = eigenmodes["eigvals"] .+ ( -J/a  + 1/tau  + J*a)

display(B[1])
display(4*A[1]- B[1]^2)
display(B[1]^2)


select = B .> 2* A

A =   A[select]
B = B[select]
theory_amin = @. pi * (B - sqrt(B^2 - 4 * A^2))/B/sqrt(-4*A^2 + 2 * B * (B - sqrt(B^2-4*A^2))) * v0^2 *2/tau*prefactor


lines!(ax, ωs[select], theory_amin, alpha=1, color="red", label="empirically exc mode - > a")

## empirical a
a =  vmean/v0

A =  sqrt.((1/tau + J * a) .* eigenmodes["eigvals"])

B = eigenmodes["eigvals"] .+ ( -J/a  + 1/tau  + J*a)

display(B[1])
display(4*A[1]- B[1]^2)
display(B[1]^2)


select = B .> 2* A

A =   A[select]
B = B[select]

theory_aemp = @. pi * (B - sqrt(B^2 - 4 * A^2))/B/sqrt(-4*A^2 + 2 * B * (B - sqrt(B^2-4*A^2))) * v0^2 *2/tau*prefactor


tag = get_tag()
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)




lines!(ax, ωs[select], theory_aemp, alpha=1, color="orange", label="empirical v/v0 - > a")

f[1,2]=Legend(f,ax)
#save(joinpath(figure_save_folder,"vprojs_projections.pdf"),f)
#xlims!(0,1)
#lines!(ax, ωs, theory_chiral, alpha=1, color = "orange")

#lines!(ax, ωs[select], theory_433, alpha=1, color = "green")
display(f)
end



    CairoMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L"\sum_{\rho}\lambda_{\rho}  (a_{\rho})^2", yscale=log10)
for i=1:10
    #lines!(ax, t, d_projs[i,:].^2)
end
lines!(ax, t, sum(d_projs.^2,dims=1 )[1,:], color="black", label=L"\sum_{\rho} (a_{\rho})^2")
lines!(ax, t, sum(eigenmodes["eigvals"] .* d_projs.^2,dims=1
 )[1,:], color="red", label=L"\sum_{\rho}\lambda_{\rho}  (a_{\rho})^2")
tag = get_tag()
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

f[1,2]=Legend(f,ax)
save(joinpath(figure_save_folder,"sum_of_mode_amplitudes_squared_times_eigenvalues.pdf"),f)
display(f)
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L"\gamma_{\nu}(v)")
vm = mean(sqrt.(vx.^2 + vy.^2),dims=1)[1,:]
vmm = mean(vm[100:end])
interval=20
maxn=100
scatter!(ax, t[500:end], vm[500:end])
for (n,eigval) in pairs(eigenmodes["eigvals"][1:interval:maxn])
    lines!(ax,t[500:end], eigval .- J*v0 ./vm[500:end] .+ J .*vm[500:end] /v0,color=1+(n-1)*interval, label="mode number $(1+(n-1)*interval)", colorrange=(1,maxn))
end

tag = get_tag()
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

f[1,2]=Legend(f,ax)
save("mode_dependent_damping_term.pdf",f)
display(f)
end


begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"\nu", ylabel=L"\gamma_{\nu}(v)")
vm = mean(sqrt.(vx.^2 + vy.^2),dims=1)[1,:]
vmm = mean(vm[4000:end])
interval=1
maxn=200
for (n,eigval) in pairs(eigenmodes["eigvals"][1:interval:maxn])

    damping =  eigval .- J*v0 ./vmm .+ J .*vmm /v0

    if damping<=0
        scatter!(ax,n,damping, color="green")

    else
        scatter!(ax,n,damping, color="red")
    end
end

tag = get_tag()
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

#f[1,2]=Legend(f,ax)
#save("mode_number_dependent_damping_term.pdf",f)
display(f)
end


begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L"\gamma_{\nu}(v)")
vm = mean(sqrt.(vx.^2 + vy.^2),dims=1)[1,:]
interval=1
maxn=41
scatter!(ax, t[500:2000], vm[500:2000]*30, label="mean v")
for (n,eigval) in pairs(eigenmodes["eigvals"][1:interval:maxn])
    lines!(ax,t[500:2000], v_projs[n,500:2000],color=1+(n-1)*interval, label="mode number $(1+(n-1)*interval)", colorrange=(1,maxn))
end

lines!(ax, t[500:2000],0.1* sin.(0.0265*t)[500:2000])

#lines!(ax, )
tag = get_tag()
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

f[1,2]=Legend(f,ax)
save("t_mode.pdf",f)
display(f)
end


CairoMakie.activate!()

min_t_ind= 500

begin
    
    f = Figure()
    ax = Axis(f[1,1], xlabel="r", ylabel= "<v(r, t)>_t/v_0")

    vpm = mean(sqrt.(vx.^2 + vy.^2)[:,min_t_ind:end],dims=2)[:,1]./v0

    rpm = mean(sqrt.(x.^2 + y.^2)[:,min_t_ind:end],dims=2)[:,1]

    scatter!(ax, rpm, vpm)

    hlines!(ax, v0/v0, color="red", label="v_0/v_0")
    hlines!(ax, vmean/v0, color="black", label="a")

    tag = get_tag()
    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
    f[1,2] = Legend(f,ax)
    ylims!(ax,0,1.1)
    display(f)
    save(joinpath(figure_save_folder,"speed_as_function_of_center_distance.pdf"),f)
end


FT_px = temporal_Fourier_transform(dt, px, min_t_ind = min_t_ind, output_not_avg=true)

CairoMakie.activate!()
begin
    
    f = Figure()
    ax = Axis(f[1,1], xlabel="r", ylabel= "ω_{px}(r)")

    rpm = mean(sqrt.(x.^2 + y.^2)[:,min_t_ind:end],dims=2)[:,1]

    scatter!(ax, rpm, FT_px["ω_max"])


    tag = get_tag()
    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
    #f[1,2] = Legend(f,ax)
    display(f)
    save(joinpath(figure_save_folder,"px_frequency_as_function_of_center_distance.pdf"),f)

end









#%% Above is already available in default analysis code
min_t_ind = 500
dt = t[2] - t[1]
FT_v_projs = temporal_Fourier_transform(dt, v_projs, min_t_ind = min_t_ind, output_not_avg=true)
FT_p_projs = temporal_Fourier_transform(dt, p_projs, min_t_ind = min_t_ind, output_not_avg=true)



CairoMakie.activate!()
GLMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"ω_{ν}", ylabel=L"FT(v_{proj})^2", yscale=log10)#, xscale=log10)
#xlims!(ax, (0,.2))
interval=2
maxn=50
tau = 1/Dr
a_min = sqrt(1 +  (eigenmodes["eigvals"][1]/(2 * J) + 1/(2 * tau *J ))^2 ) - (eigenmodes["eigvals"][1]/(2 * J) + 1/(2 * tau *J))


damping = eigenmodes["eigvals"][1] + J * a_min * (1 - 1/a_min^2 )+1/tau
display(damping)

display("a_min = $a_min")
#scatter!(ax, t[500:2000], vm[500:2000]*30)
for (n,eigval) in pairs(eigenmodes["eigvals"][1:interval:maxn])
    scatter!(ax, FT_v_projs["ω"], FT_v_projs["Xf2"][n,:],color=1+(n-1)*interval, colorrange=(1,maxn),label="mode $(1+(n-1)*interval)",alpha=0.1)



    tau = 1/Dr

    
    #display(a_min)
    
    #a = a_min
   # display(a)
    
    ω = FT_v_projs["ω"][1:end]

    theory_ABP = @.  v0^2 * 2*tau * 2* pi *  ω^2 / ((1+(tau*ω)^2)* (eigval^2 +  ω^2))

    a =  vmean/v0


    A =  sqrt.((1/tau + J * a) .* eigval)

    B = eigval .+ ( -J/a  + 1/tau  + J*a)

    theory_J = @. 2*pi * v0^2 * 2/tau * ω^2/( (A^2 - ω^2)^2. + ω^2 * B^2)

    # theory =ω.^2 .* ( 1 .+  2 * J/a * eigval ./ (ω.^2 .+ eigval^2) .+ 0* J^2 * a^2 * (eigval^2 .+ ω.^2 .* (1 - 1/a^2)^2) ./ (ω.^4 + ω.^2 .* eigval^2)) .* 1 ./(eigval^2 .+ ω.^2) .*2*pi*1/Dr * v0^2 ./(1 .+ (ω ./ Dr).^2) *  (length(t)-min_t_ind)/(2*pi)
    
    # theory_J0 = ω.^2 .* ( 1) .* 1 ./(eigval^2 .+ ω.^2) .*2*pi*1/Dr * v0^2 ./(1 .+ (ω ./ Dr).^2) *  (length(t)-min_t_ind)/(2*pi)

    # subterm = (eigval + J * a * (1 - 1/a^2))

    # subterm = (eigval )
    # denominator = (ω.^2 .-  J * a * eigval).^2 .+ ω.^2 .* subterm^2


    #Np comes from the assumption that the noise vectors are to replaced by sum of iid gaussian noise vectors

    # 1e-5 to offset the plot

    #elastic_offset =0.0
    #theory_433 =   ω.^2 .* v0^2 * 2/tau * 2 * pi^2 ./ ( ( (1/tau + J *( a + elastic_offset) )*eigval .- ω.^2).^2 .+ ω.^2 .* (eigval - J/a + 1/tau + J *(a+elastic_offset) )^2)*  ((length(t)-min_t_ind)/2/pi)

    
    #theory_433_no_J =  ω.^2 .* v0^2 * 2/tau * 2 * pi ./ ( ( (1/tau + 0* J * a)*eigval .- ω.^2).^2 .+ ω.^2 .* (eigval - 0* J/a + 1/tau + 0*J * a)^2)*  (length(t)-min_t_ind)/(2*pi)


    #hlines!(ax, (eigval_s - J/a + 1/tau + J * a)^2,colorrange=(1,maxn) ,color=1+(n-1)*interval)

    #theory_v2 = ω.^2 .* ( 1 .+1* 2* J * a * ( eigval .* ω.^2  .- ω.^2 .* (1 - 1/a^2) .* subterm .- J* a * eigval^2  )./ denominator .+1*  J^2 * a^2*(eigval^2  .+ ω.^2 .*  (1 - 1/a^2)^2 )./denominator  ).* 1 ./(eigval^2 .+ ω.^2) .*2*pi*1/Dr * v0^2 ./(1 .+ (ω ./ Dr).^2) *  (length(t)-min_t_ind)/(2*pi)


    
    #lines!(ω, theory,color=1+(n-1)*interval, colorrange=(1,maxn),label="mode $(1+(n-1)*interval)")

    #lines!(ax,ω, theory_433 ,color=1+(n-1)*interval, colorrange=(1,maxn),label="mode $(1+(n-1)*interval)", linestyle=:dash)


    #vlines!(ax, sqrt((1/tau + J * a)*eigval),color=1+(n-1)*interval, colorrange=(1,maxn),label="mode $(1+(n-1)*interval)", linestyle=:dash)
    lines!(ω, theory_ABP*  (length(t)-min_t_ind)/(2*pi) ,color=1+(n-1)*interval, colorrange=(1,maxn),label="mode $(1+(n-1)*interval)", linestyle=:dash,alpha=0.2)
    lines!(ω, theory_J*  (length(t)-min_t_ind)/(2*pi) ,color=1+(n-1)*interval, colorrange=(1,maxn),label="mode $(1+(n-1)*interval)", linestyle=:solid, alpha=0.2)
end

xlims!(ax, (0.001,.2))
tag = get_tag()
#vlines!(ax, FT_px["ω_max"][1], label="ω_{px}", color="black")
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

f[1,2]=Legend(f,ax)
#save(joinpath(figure_save_folder,"Np_mode_v_proj_frequency_prediction.pdf"),f)
display(f)
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"ω_{ν}", ylabel=L"FT(v_{proj})^2", yscale=log10)
xlims!(ax, (0,.2))
interval=1
maxn=50
vmm=mean(vm[500:end])
#scatter!(ax, t[500:2000], vm[500:2000]*30)
for (n,eigval) in pairs(eigenmodes["eigvals"][1:maxn])

    if n%1===0 && n < maxn
        scatter!(ax, FT_v_projs["ω"], FT_v_projs["Xf2"][n,:],color=1+(n-1)*interval, colorrange=(1,maxn),label="mode $(1+(n-1)*interval)")

        γ = eigval - J*v0/vmm + J*vmm/v0

        Ω2  = J *vmm/v0 * eigval

        prediction = sqrt( J * ( sqrt( 1+(eigval/(2*J))^2 ) - eigval/(2*J) ) * eigval )


        #vlines!(ax, 2*abs(sqrt( Complex((γ/2)^2 -Ω2 ) ) ),color=1+(n-1)*interval, colorrange=(1,maxn))

        vlines!(ax, prediction,color=1+(n-1)*interval, colorrange=(1,maxn))
    end
    #vlines!(ax, -eigval^2 + J*v0/vmm),color="red")#1+(n-1)*interval, colorrange=(1,maxn))
    #vlines!(ax, 2 * sqrt(Ω2) ,color=1+(n-1)*interval, colorrange=(1,maxn))
end

#vlines!(ax,  sqrt( J*eigenmodes["eigvals"][40]*( sqrt(1+(eigenmodes["eigvals"][40]/(2*J))^2 ) - eigenmodes["eigvals"][40]/(2*J)) ),color="red")

tag = get_tag()
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

f[1,2]=Legend(f,ax)
#save("Np_mode_v_proj_frequency_prediction.pdf",f)
display(f)
end







FT = temporal_Fourier_transform(t[2]-t[1],px,min_t_ind=500)

begin
f = Figure()#
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{p_x(t)\}(\omega)|^2");
scatter!(ax, FT["ω"], FT["pavg_X2"])
scatter!(ax,FT["ω_max"] , FT["max_X2"],label="max")


vlines!(ax, J*sqrt(1 - vmean^2/v0^2),color="red")
Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)
xlims!(ax, 0,0.2)

#vlines!(ax, sqrt(J*k*( sqrt(1+(k/(2*J))^2 ) - k/(2*J))), label="theory", color="green")
f[1,2]=Legend(f,ax)
#ylims!(ax, low=1e-7,high= 1e2)
#save("single_particle_small_noise_unit_alignment.pdf",f)
display(f)
end
begin

Nint = size(vx)[1]
it = 1500

B = zeros( 2* Nint, 2*Nint)
C = zeros( 2* Nint, 2*Nint)


for i in 1:2*Nint


    B[i,i] = 1/vit[i]
    C[i,i] = vit[i]

end

speedit = sqrt.( vx[:,it].^2 +  vy[:,it].^2)
vit =collect(Iterators.flatten(zip(speedit, speedit)))

begin
f = Figure()
ax = Axis(f[1,1]);

hist!(ax, speedit/v0, normalization = :pdf)
xlims!(ax, 0,2)
display(f)
end
end



dRDdVt = transpose((D * collect(Iterators.flatten(zip(vx[:,it], vy[:,it]))) )) * collect(Iterators.flatten(zip(dx[:,it], dy[:,it])))
M = D - J*v0*B + J/v0 * C + J/v0 * B * dRDdVt

begin
f = Figure()#
ax = Axis(f[1,1]);
heatmap!(ax, D - J*v0*B + J/v0 * C, colorrange=(-1,1))

display(f)
    
end
begin
f = Figure()#
ax = Axis(f[1,1]);
heatmap!(ax, D)

display(f)
    
end
begin
f = Figure()#
ax = Axis(f[1,1]);
heatmap!(ax, B, colorrange=(-1,1))

display(f)
    
end


begin
it=4000
speedit = sqrt.( vx[:,it].^2 +  vy[:,it].^2)
f = Figure()#
ax = Axis(f[1,1]);
scatter!(ax, x[:,it], y[:,it], color = speedit)

display(f)
    

    
end


begin
vals, inds = findmax(FT_v_projs["Xf2"], dims=2)
max_ωs = FT_v_projs["ω"][getindex.(inds,2)][:,1]

f = Figure()
ax = Axis(f[1,1], xlabel=L"ω_{ν}", ylabel=L" ω_{max}")
scatter!(ax, sqrt.(eigenmodes["eigvals"]),max_ωs)

display(f)
end


begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"ω_{ν}", ylabel=L"\omega")
modenumbers = range(1,size(v_projs)[1])
heatmap!(ax, sqrt.(eigenmodes["eigvals"]),FT_v_projs["ω"], log10.(FT_v_projs["Xf2"]))
display(f)
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"ω_{ν}", ylabel=L"\omega")
heatmap!(ax, sqrt.(eigenmodes["eigvals"]),FT_p_projs["ω"], log10.(FT_p_projs["Xf2"]))

display(f)
end


vrms = sqrt( mean(vx.^2 + vy.^2) )


begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L"p_proj")
scatter!(ax, t, p_projs[1,:])
lines!(ax, t, sin.(t * 0.0266))
scatter!(ax, t, mean(px, dims=1)[1,:])
display(f)
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L"p_proj")

lines!(ax, t, sum(v_projs.^2, dims=1)[1,:])
display(f)
end

θp1 = angle.( px[1,:] + im * py[1,:])


using GLMakie
GLMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L" θ_p")
#lines!(ax, t, px[1,:])
scatterlines!(ax, t, θp1)
display(f)
end

d= diff(θp1)

θp1c =  [0.]
for(i, di) in pairs(d) 

    if -pi<di<pi
        push!(θp1c, θp1c[i] + di)

    elseif di>pi
        push!(θp1c, θp1c[i] + di-2*pi)

    elseif di<-pi
        push!(θp1c, θp1c[i] + di+2*pi)
    end
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L" θ_p")
#lines!(ax, t, px[1,:])
scatterlines!(ax, t, θp1c)
scatterlines!(ax, t, θp1)
display(f)
end



θpc = unwrap(angle.( px + im * py))

θvc = unwrap(angle.( vx + im * vy))




begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L" θ_p")
#lines!(ax, t, px[1,:])

for p in 1:size(θpc)[1]
    if p==1000
    #scatterlines!(ax, t, θvc[p,:])
    #scatterlines!(ax, t, θpc[p,:])
    scatterlines!(ax, t, θvc[p,:] -θpc[p,:])
    end
end
tag = get_tag()
#f[1,2]=Legend(f,ax)
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
display(f)

end



using CairoMakie
CairoMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L" θ_p")
#lines!(ax, t, px[1,:])

for p in 1:size(θpc)[1]
    if p==1
    scatterlines!(ax, t, θpc[p,:], label="unwrapped θ_p")

    scatter!(ax, t, cos.( θpc[p,:]), markersize=10, color="orange",label="cos θ_p")
    lines!(ax, t, px[p,:], color="green",label="px")

    end
end
tag = get_tag()
f[1,2]=Legend(f,ax)
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
display(f)
save("unwrap_angle.pdf",f)
end



begin

msd = zeros(length(t))
for (i, ti) in pairs(t)

    msd[i] = mean(  (θpc[1,i] -θpc[1,1]) .^2 )
end
end
begin
using GLMakie 
GLMakie.activate!()
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L" 1", yscale=log10)
lines!(ax, t, msd)
display(f)
end


begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L" v")
lines!(ax, t, sqrt.(vx.^2 + vy.^2)[1,:])

lines!(ax, t, vx[1,:])
display(f)
end





begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=L" 1")
lines!(ax, t, px[1,:])
lines!(ax, t, vx[1,:]/v0)
display(f)
end


using CairoMakie
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"t", ylabel=" x-component averaged over particles")
lines!(ax, t, mean(px,dims=1)[1,:], label="px")
lines!(ax, t,mean(vx,dims=1)[1,:]/v0, label="vx/v0")
tag = get_tag()
f[1,2]=Legend(f,ax)
Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
display(f)
save("x-components.pdf",f)
end



begin
    f = Figure()
    ax = Axis(f[1,1], xlabel=L"px", ylabel="vx")

    for (i, ti) in pairs(t)
        scatter!(ax, px[:,i], vx[:,i], color=i, colorrange=(0, t[end]))
    end


    tag = get_tag()
    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
    display(f)  

end

begin
    f = Figure()
    ax = Axis(f[1,1], xlabel=L"px", ylabel="py")

    for (i, ti) in pairs(t)
        scatter!(ax, px[:,i], py[:,i], color=i, colorrange=(0, t[end]))
    end


    tag = get_tag()
    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
    display(f)  

end

begin
    f = Figure()
    ax = Axis(f[1,1], xlabel=L"vx", ylabel="vy")

    for (i, ti) in pairs(t)
        scatter!(ax, vx[:,i], vy[:,i], color=i, colorrange=(0, t[end]))
    end


    tag = get_tag()
    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
    display(f)  

end



begin
    using CairoMakie
    CairoMakie.activate!() 
    f = Figure(size=(5000,5000))
    ax = Axis(f[1,1], xlabel=L"x", ylabel="y", aspect=1)

    for (i, ti) in pairs(t)

        if i%10==0
        scatter!(ax, x[:,i], y[:,i], color=i, colorrange=(0, t[end]), colormap=:viridis)
        end
    end


    tag = get_tag()
    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
    display(f)  

end



