
include(joinpath("..","src","Engine.jl"))
include("AnalysisPipeline.jl")
include("AnalysisFunctions.jl")


#base_folder = joinpath(homedir(),"sa","survey","hex_disordered","phi_1","Nlin_4","vary_J_Dr")

base_folder = joinpath("/Volumes","T7_Shield","sa","single","noise","not_all_1")

#base_folder = joinpath(homedir(),"sa","survey","hex_disordered","phi_1","Nlin_4","vary_J_Dr")
raw_data_base_folder = joinpath(base_folder, "simdata")

raw_data_file_path = joinpath(raw_data_base_folder,"raw_data.jld2")

raw_data_file = jldopen(raw_data_file_path, "r")

frames = raw_data_file["frames"]
system = raw_data_file["system"]

#frames_support = jldopen(joinpath(raw_data_base_folder,"J_0.1","Dr_0.01","seed_23","ra_raw_data.jld2"),"r")["frames"]

t =  raw_data_file["integration_info"]["save_tax"]
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

@views for i in 1:Nt
    x[:,i] = extract_frame_data_for_type("x", 1, frames[string(i)])
    y[:,i] = extract_frame_data_for_type("y", 1, frames[string(i)])



    vx[:,i] = extract_frame_data_for_type("vx", 1, frames[string(i)])
    vy[:,i] = extract_frame_data_for_type("vy", 1, frames[string(i)])

    px[:,i] = extract_frame_data_for_type("px", 1, frames[string(i)])
    py[:,i] = extract_frame_data_for_type("py", 1, frames[string(i)])

    qx[:,i] = extract_frame_data_for_type("qx", 1, frames[string(i)])
    qy[:,i] = extract_frame_data_for_type("qy", 1, frames[string(i)])

end

FT = temporal_Fourier_transform(t[2]-t[1],px,min_t_ind=100)




using CairoMakie
CairoMakie.activate!()
begin
f = Figure()
ax = Axis(f[1,1],yscale=log10, xlabel="ω",ylabel=L"|  \mathcal{F}\{p_x(t)\}(\omega)|^2");
scatter!(ax, FT["ω"], FT["pavg_X2"])
scatter!(ax,FT["ω_max"] , FT["max_X2"],label="max")



Label(f[2,1],"Dr = $Dr, J = $J, v_0 = $v0, k=$k, ", tellwidth=false, halign=:left, word_wrap = true)


vlines!(ax, sqrt(J*k*( sqrt(1+(k/(2*J))^2 ) - k/(2*J))), label="theory", color="green")
f[1,2]=Legend(f,ax)
#ylims!(ax, low=1e-7,high= 1e2)
save("single_particle_small_noise_unit_alignment.pdf",f)
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
save("single_particle_small_noise_unit_alignment_zoomed.pdf",f)
display(f)
end