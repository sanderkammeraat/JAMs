include(joinpath("..","src","Engine.jl"))
include("AnalysisPipeline.jl")
include("AnalysisFunctions.jl")


#base_folder = joinpath(homedir(),"sa","survey","hex_disordered","phi_1","Nlin_4","vary_J_Dr")

base_folder = joinpath("/Volumes","T7_Shield","sa","survey","hex_disordered", "phi_1", "Nlin_20", "vary_J_Dr")

#base_folder = joinpath(homedir(),"sa","survey","hex_disordered","phi_1","Nlin_4","vary_J_Dr")
raw_data_base_folder = joinpath(base_folder, "simdata")

#Make tree to navigate simulation data folder structure
tree = construct_folder_tree_param_param_seed(raw_data_base_folder)

raw_data_file_path = joinpath(raw_data_base_folder,"J_0.1","Dr_0.01","seed_23","sa_raw_data.jld2")

raw_data_file = jldopen(raw_data_file_path, "r")

frames = raw_data_file["frames"]
system = raw_data_file["system"]

frames_support = jldopen(joinpath(raw_data_base_folder,"J_0.1","Dr_0.01","seed_23","ra_raw_data.jld2"),"r")["frames"]

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
    x[:,i] = extract_frame_data_for_type("x", 1, frames[string(i)])
    y[:,i] = extract_frame_data_for_type("y", 1, frames[string(i)])



    vx[:,i] = extract_frame_data_for_type("vx", 1, frames[string(i)])
    vy[:,i] = extract_frame_data_for_type("vy", 1, frames[string(i)])

    px[:,i] = extract_frame_data_for_type("px", 1, frames[string(i)])
    py[:,i] = extract_frame_data_for_type("py", 1, frames[string(i)])

    qx[:,i] = extract_frame_data_for_type("qx", 1, frames[string(i)])
    qy[:,i] = extract_frame_data_for_type("qy", 1, frames[string(i)])

end


x0 = frames_support[string(length(frames_support))]["x"]
y0 = frames_support[string(length(frames_support))]["y"]

x0int = extract_frame_data_for_type("x",1,frames_support[string(length(frames_support))])
y0int = extract_frame_data_for_type("y",1,frames_support[string(length(frames_support))])
#With boundary
D_wb = construct_D(x0, y0, k, R, type)

#Boundary particles are stored at the end
interior_indmax = sum(type.==1)
D = Symmetric(D_wb[1:2*interior_indmax,1:2*interior_indmax])


modes = diagonalize_D(D)



v_projs = project_on_eigvecs(modes["eigvecs"], vx,vy)
p_projs = project_on_eigvecs(modes["eigvecs"], px,py)

#%% Above is already available in default analysis code

dt = t[2] - t[1]
FT_v_projs = temporal_Fourier_transform(dt, v_projs, min_t_ind = 500, output_not_avg=true)
FT_p_projs = temporal_Fourier_transform(dt, p_projs, min_t_ind = 500, output_not_avg=true)
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"ω_{ν}", ylabel=L"\omega")
modenumbers = range(1,size(v_projs)[1])
heatmap!(ax, sqrt.(modes["eigvals"]),FT_v_projs["ω"], log10.(FT_v_projs["Xf2"]))
display(f)
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"ω_{ν}", ylabel=L"\omega")
heatmap!(ax, sqrt.(modes["eigvals"]),FT_p_projs["ω"], log10.(FT_p_projs["Xf2"]))

display(f)
end

vals, inds = findmax(FT_v_projs["Xf2"], dims=2)
max_ωs = FT_v_projs["ω"][getindex.(inds,2)][:,1]

vrms = sqrt( mean(vx.^2 + vy.^2) )
begin
f = Figure()
ax = Axis(f[1,1], xlabel=L"ω_{ν}", ylabel=L" ω_{max}")
scatter!(ax, sqrt.(modes["eigvals"]),max_ωs)

display(f)
end

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



