




include(joinpath("..","src","Engine.jl"))
include("AnalysisPipeline.jl")
include("AnalysisFunctions.jl")


tree = construct_folder_tree_param_param_seed(joinpath(homedir(),"sa","survey","hex_disordered", "phi_1","Nlin_4","vary_J_Dr","simdata"))


#tree = construct_folder_tree_param_param_seed(joinpath(homedir(),"mounting","data1_kammeraat","sa","survey","hex_disordered", "phi_1","Nlin_20","vary_J_Dr","simdata"))

data = jldopen( joinpath(tree["J_1.0"]["Dr_0.01"]["seed_53"], "sa_raw_data.jld2") )

#data = jldopen(joinpath(homedir(),"mounting","data1_kammeraat","sa", "phi_1","Nlin_20","vary_J_Dr","simdata","J_0.1","Dr_0.01","seed_7","raw_data.jld2"))

#Ignore boundary particles
frames = data["frames"]
t = data["integration_info"]["save_tax"]

Nt = length(t)
Np = length(extract_frame_data_for_type("id",1,frames["1"]))
x = zeros(Np, Nt)
y = zeros(Np, Nt)
vx = zeros(Np, Nt)
vy = zeros(Np, Nt)

print(data["system"])

px = zeros(Np, Nt)
py = zeros(Np, Nt)

@views for i in 1:Nt
    x[:,i] = extract_frame_data_for_type("x", 1, frames[string(i)])
    y[:,i] = extract_frame_data_for_type("y", 1, frames[string(i)])

    vx[:,i] = extract_frame_data_for_type("vx", 1, frames[string(i)])
    vy[:,i] = extract_frame_data_for_type("vy", 1, frames[string(i)])

    px[:,i] = extract_frame_data_for_type("px", 1, frames[string(i)])
    py[:,i] = extract_frame_data_for_type("py", 1, frames[string(i)])

end

dx = x[:,2:end] .- x[:,1]

dy = y[:,2:end] .- y[:,1]

dr = sqrt.(dx.^2 + dy.^2)

using GLMakie
GLMakie.activate!()

begin
f = Figure();
fidmax = 1000
pid = 10
ax = Axis(f[1,1], aspect=DataAspect())

scatterlines!(ax, dx[pid,1:fidmax], dy[pid,1:fidmax], color = t[1:fidmax], colormap=:rainbow)
display(f)          
end

begin
    f = Figure();
    fidmax = 1000
    pid = 10
    ax = Axis(f[1,1])
    
    scatterlines!(ax,t[1:fidmax], dx[pid,1:fidmax])
    display(f)          
end




begin
    f = Figure();
    fidmax = 1000
    pid = 10
    ax = Axis(f[1,1])
    
    scatterlines!(ax,t[2:end], (dx./dr)[pid,:], color = t[2:end], colormap=:rainbow)
    display(f)

end


function p_correlation(binsize, maxbin_center)


    rbin_edges = prepend!(collect(range(start=0, step=binsize, stop=maxbin_center)),[0])

    rbin_edges2 = rbin_edges.^2

    rbin_centers = (rbin_edges[2:end] + rbin_edges[1:end-1])/2

    Nbin = length(rbin_centers)

    C = ones(size(px)[2], Nbin)*NaN

    counts = zeros(size(px)[2], Nbin)

    for i in eachindex(t)
        for p1 in 1:size(px)[1]

            for p2 in p1:size(px)[1]

                Δr2 = (x[p1, i]- x[p2,i])^2 +   (y[p1, i]- y[p2,i])^2

                for bin in eachindex(rbin_centers)

                    if (Δr2<= rbin_edges2[bin+1] && Δr2> rbin_edges2[bin]) || (p1==p2 && bin==1)

                        if isnan(C[i, bin])
                            C[i, bin]=0
                            counts[i, bin]=0
                        end

                        C[i, bin] += px[p1, i]* px[p2, i] + py[p1, i] * py[p2, i]
                        counts[i,bin] += 1
                    end

                end

            end
        end
    end

    C.=C./counts
    return rbin_centers, C

end

PCOR = spatial_p_correlation(2.5, 10, x,y,px, py)


@time spatial_p_correlation(2.5, 10, x,y,px, py)

begin
    f = Figure();
    ax = Axis(f[1,1])
    
    scatterlines!(ax,PCOR["rbc"], mean(PCOR["C"][500:end,:], dims=1)[1,:] )
    display(f)

end


PSTCOR = spatiotemporal_p_correlation(2.5, 10, x[:,1],y[:,1],px, py, min_t_ind=500)
begin
    f = Figure();
    ax = Axis(f[1,1], yreversed=true)
    t_max_ind = 5000
    
    heatmap!(ax,  PSTCOR["rbe"],t[1:t_max_ind],   transpose(PSTCOR["C"][1:t_max_ind,:]), colormap=:seismic, colorrange=(-1,1))
    vlines!(ax, PSTCOR["rbe"], color=:black)
    display(f)

end

@profview spatiotemporal_p_correlation(2.5, 10, x[:,1],y[:,1],px, py)


begin
    f = Figure();
    ax = Axis(f[1,1])
    t_max_ind = 5000
    scatterlines!(ax,  PSTCOR["rbc"],PSTCOR["C"][1,:])
    scatterlines!(ax,  PSTCOR["rbc"],PSTCOR["C"][end-497,:])
    display(f)

end

xmax = maximum(x)
ymax = maximum(y)

rmax = sqrt(xmax^2 + ymax^2)

using JLD2

analysis = jldopen("/Users/kammeraat/sa/phi_1/Nlin_4/vary_J_Dr/analysis_v2/J_0.01/Dr_0.0/seed_10.jld2")

analysis = jldopen("/Users/kammeraat/mounting/data1_kammeraat/sa/phi_1/Nlin_20/vary_J_Dr/analysis_v2/J_1.0/Dr_1.0/seed_14.jld2")

SPTE_p=analysis["SPTE_p"]
t = analysis["t"]
begin
    f = Figure();
    ax = Axis(f[1,1], yreversed=true)
    t_max_ind = 5000
    
    heatmap!(ax,  SPTE_p["rbe"],t[1:t_max_ind],   transpose(SPTE_p["C"][1:t_max_ind,:]), colormap=:seismic, colorrange=(-1,1))
    vlines!(ax, SPTE_p["rbe"], color=:black)
    display(f)

end


#%%


using FFTW

tss = 500

FT =rfft(dx[:,tss:end],2)
FTfreq = rfftfreq( size(dx[:,tss:end])[2], 1/ (t[2]-t[1])).*2*pi


begin
f = Figure();

ax = Axis(f[1,1], yscale=log10);
for i in 1:size(FT)[1]

scatter!(ax, FTfreq, abs.(FT[i,:]).^2, color=i, colorrange=(1,size(FT)[1] ))

end
display(f)

end

FT_particle_mean = mean(abs.(FT).^2, dims=1)[1,:]

maxval, maxind = findmax(FT_particle_mean)

main_freq = FTfreq[maxind]

begin
    f = Figure();
    
    ax = Axis(f[1,1], yscale=log10, xlabel=L"\omega", ylabel=L"|  \mathcal{F}\{\delta x(t)\}(\omega)|^2");
    
    
    scatter!(ax, FTfreq, FT_particle_mean)

    scatter!(ax,main_freq , maxval)
    
    display(f)
    
end

2*pi/main_freq

using Peaks

begin
f = Figure();
    
ax = Axis(f[1,1], yscale=log10, xlabel=L"\omega", ylabel=L"|  \mathcal{F}\{\delta x(t)\}(\omega)|^2");
inds, vals = findmaxima(FT_particle_mean,200)
scatter!(ax, FTfreq, FT_particle_mean)
scatter!(ax,FTfreq[inds] , vals)

display(f)
end

FTfreq[inds][2]/FTfreq[inds][1]





FT = temporal_Fourier_transform(t[2]-t[1], dx, min_t_ind=500)


begin
    f = Figure();
        
    ax = Axis(f[1,1], yscale=log10, xlabel=L"\omega", ylabel=L"|  \mathcal{F}\{\delta x(t)\}(\omega)|^2");
    inds, vals = findmaxima(FT_particle_mean,200)
    band!(ax, FT["ω"], FT["pavg_X2"] .-FT["pste_X2"] ,FT["pavg_X2"] .+FT["pste_X2"], color=:red)
    scatter!(ax, FT["ω"], FT["pavg_X2"])
    scatter!(ax,FT["ω_max"] , FT["max_X2"])


    
    display(f)
end




close(data)