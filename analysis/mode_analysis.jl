include(joinpath("..","src","Engine.jl"))

include(joinpath("LinearizeSystem.jl"))

using CairoMakie #vector graphics
begin #load and check frames

base_folder = joinpath(homedir(),"sa","production","phi_1","tstop_2e3")

rx_folder = joinpath(base_folder, "simdata","relaxation_step")

sa_folder = joinpath(base_folder, "simdata", "self_alignment_step")

ra_folder = joinpath(base_folder, "simdata", "relax_again_step")
rx_file = jldopen( joinpath(rx_folder, "raw_data.jld2"), "r")
sa_file = jldopen( joinpath(sa_folder, "raw_data.jld2"), "r")
ra_file = jldopen( joinpath(ra_folder, "raw_data.jld2"), "r")

rx_final_frame = rx_file["frames"][string(length(rx_file["frames"]))]

ra_final_frame = ra_file["frames"][string(length(ra_file["frames"]))]

#Is final relaxation step frame x,v same as initial self-alignment step frame x, y
rx_xf = rx_final_frame["x"]

rx_tf = rx_final_frame["t"]


sa_first_frame = sa_file["frames"]["1"]

sa_xi = sa_first_frame["x"]

sa_ti = sa_first_frame["t"]

f,ax,=plot(sa_xi)
plot!(ax,rx_xf, marker=:diamond, markersize=6)
@assert sa_xi == rx_xf


sa_yi = sa_first_frame["y"]

#So let's take

x0_all = sa_xi
y0_all = sa_yi
R_all = sa_first_frame["R"]
type_all= sa_first_frame["type"]


f = Figure()
ax=Axis(f[1,1], xlabel = "x", ylabel="y", aspect=1)
scatter!(ax,x0_all,y0_all, markerspace=:data, markersize=R_all, color = type_all)
display(f)

end

# for (i,reference_frame) in pairs([sa_first_frame, ra_final_frame])

#     if i==1
#         name = "using_sa_first_frame_as_reference"

#     else
#         name = "using_ra_final_frame_as_reference"
#     end
#     figure_base_folder=mkpath(joinpath(base_folder,name))


reference_frame = ra_final_frame

x0 =reference_frame["x"]


y0 =reference_frame["y"]

R= reference_frame["R"]

type= reference_frame["type"]
id=reference_frame["id"] 

f = Figure();
ax=Axis(f[1,1], xlabel = "x", ylabel="y", aspect=1)
scatter!(ax,x0,y0, markerspace=:data, markersize=R, color = type)


k = ra_file["system"]["forces"]["pair_forces"]["soft_disk_force"]["karray"]


D_wb = construct_D(x0,y0,k,R,type)

interior_indmax = sum(type.==1)

D = Symmetric(D_wb[1:2*interior_indmax,1:2*interior_indmax])


f = Figure();
ax=Axis(f[1,1], xlabel = "2j", ylabel="2i", title="D entries", yreversed=true,aspect = DataAspect())
hm =heatmap!(ax, D, colormap=:viridis)
Colorbar(f[:, end+1], hm)
display(f);
save(joinpath(figure_base_folder,"D.pdf"),f)







eigenfact = eigen(D)

eigenfact.vectors

f = Figure();
ax=Axis(f[1,1], xlabel = "m", ylabel="λ_m", aspect=1, title="Eigenvalues")
plot!(ax,eigenfact.values, markersize=1)

display(f)
save(joinpath(figure_base_folder,"eigenvalues.pdf"),f)

f = Figure(size=(4000,4000));
mm=1

for i in 1:4
    for j in 1:4
    print(mm)

    m = mm
    

    eigen_m_value = eigenfact.values[m]
    eigen_m_vector =  reshape(eigenfact.vectors[:,m], (2, round(Int64,length(eigenfact.vectors[:,m])/2)) )
    eigen_m_x = eigen_m_vector[1,:]
    eigen_m_y = eigen_m_vector[2,:]


    ax=Axis(f[i,j], xlabel = "x", ylabel="y", aspect=1, title=string(eigen_m_value))
    scatter!(ax,x0,y0, markerspace=:data, markersize=2 .*R.^2, color = id,marker = Circle,alpha=0.1)

    arrows!(ax, x0, y0, eigen_m_x, eigen_m_y, lengthscale = 100)
    mm+=1
    end

end
display(f)
save(joinpath(figure_base_folder,"modevis.pdf"),f)


f = Figure();
ax=Axis(f[1,1], ylabel = "pdf(ω)", xlabel="ω = sqrt(λ_m)", aspect=1, title="Density of states")
hist!(ax,sqrt.(eigenfact.values[eigenfact.values .>=0]), bins=50, normalization=:pdf);
#plot!(ax,sqrt.(eigenfact.values[eigenfact.values .>=0]),sqrt.(eigenfact.values[eigenfact.values .>=0]).^2/100);
display(f);
save(joinpath(figure_base_folder,"density_of_states.pdf"),f)



#Now do projections over time

Nframes = length(sa_file["frames"])
Neigenvectors =  length(eigenfact.values)
pos_projections = zeros(Nframes,Neigenvectors )
vel_projections= zeros(Nframes,Neigenvectors )

x0interior =  extract_data_for_type("x", 1, reference_frame)
y0interior =  extract_data_for_type("y", 1, reference_frame)
vrms = zeros(Nframes)
@views for i in 1:Nframes

    xi = extract_data_for_type("x", 1, sa_file["frames"][string(i)])
    yi = extract_data_for_type("y", 1, sa_file["frames"][string(i)])

    vxi = extract_data_for_type("vx", 1, sa_file["frames"][string(i)])
    vyi = extract_data_for_type("vy", 1, sa_file["frames"][string(i)])

    #Interweave and calculate normalized displacement
    displacement = collect(Iterators.flatten(zip(xi .- x0interior, yi .- y0interior)))

    #Interweave and calculate normalized velocity
    velocity = collect(Iterators.flatten(zip(vxi, vyi)))

    for j in 1:Neigenvectors
        pos_projections[i,j] = sum( displacement .* eigenfact.vectors[:,j] )
        vel_projections[i,j] = sum( velocity .* eigenfact.vectors[:,j] )
    end

    vrms[i] = sqrt(mean(vxi.^2 + vyi.^2))

end

f = Figure();
ax=Axis(f[1,1], ylabel = "<λ_m|dr(t)>", xlabel="t", title="displacement projections")
save_tax = sa_file["integration_info"]["save_tax"]
Nplot = 200
#ylims!(-0.1,0.1)
for i in  1:Nplot
    lines!(ax, save_tax, pos_projections[:,i], color=eigenfact.values[i], colorrange=(0,eigenfact.values[Nplot]), colormap=:viridis)


end
display(f)
save(joinpath(figure_base_folder,"displacement_projections.pdf"),f)

f = Figure();
ax=Axis(f[1,1], ylabel = "<λ_m|v(t)>", xlabel="t", title="velocity projections")
save_tax = sa_file["integration_info"]["save_tax"]
Nplot = 200
#ylims!(-0.1,0.1)
for i in  1:Nplot
    lines!(ax, save_tax, vel_projections[:,i], color=eigenfact.values[i], colorrange=(0,eigenfact.values[Nplot]), colormap=:viridis)


end
display(f)
save(joinpath(figure_base_folder,"velocity_projections.pdf"),f)

#end # end of loop over possible reference state

begin #Plotting angular dynamic
Nframes = length(sa_file["frames"])
figure_folder_angular_dynamics = mkpath(joinpath(homedir(), "Obsidian","Sync","figure_warehouse","angular_dynamics" ))
Np = 1000
vp_angles = zeros(Nframes, Np)

p_angles = zeros(Nframes, Np)

v_angles = zeros(Nframes, Np)

@views for i in 1:Nframes


    vxi = extract_data_for_type("vx", 1, sa_file["frames"][string(i)])
    vyi = extract_data_for_type("vy", 1, sa_file["frames"][string(i)])

    pxi = extract_data_for_type("px", 1, sa_file["frames"][string(i)])
    pyi = extract_data_for_type("py", 1, sa_file["frames"][string(i)])

    v_angles[i,:] = angle.(vxi .+1im .*vyi)
    p_angles[i,:] = angle.(pxi .+1im .*pyi)

    vp_angles[i,:]=angle.( exp.( 1im * (angle.(vxi .+1im .*vyi) - angle.(pxi .+ 1im .*pyi)  )  )  )

end

f = Figure();
ax=Axis(f[1,1], ylabel = "angle(v,p)>", xlabel="t", title="vp angles")
save_tax = sa_file["integration_info"]["save_tax"]
#ylims!(-0.1,0.1)
for i in  1:Np
    lines!(ax, save_tax, vp_angles[:,i],colormap=:viridis)

    #lines!(ax, save_tax, p_angles[:,i],colormap=:viridis)

    #lines!(ax, save_tax, vp_angles[:,i],colormap=:viridis)

    
end

pm = mean(p_angles, dims=2)[:,1]
pstd = std(p_angles, dims=2)[:,1]


vm = mean(v_angles, dims=2)[:,1]
vstd =  std(v_angles, dims=2)[:,1]


vpm =  mean(vp_angles, dims=2)[:,1]
vpstd =  std(vp_angles, dims=2)[:,1]

f = Figure(size=(2000,500));
ax=Axis(f[1,1], ylabel = "mean angles", xlabel="t", title="mean angles")

save_tax = sa_file["integration_info"]["save_tax"]
p=plot!(ax, save_tax, pm)
band!(ax, save_tax, pm-pstd, pm+pstd, alpha=0.5)


v=plot!(ax, save_tax,vm ,colormap=:viridis)
band!(ax, save_tax, vm-vstd, vm+vstd, alpha=0.5)

vp=plot!(ax, save_tax,vpm,colormap=:viridis)
band!(ax, save_tax, vpm-vpstd, vpm+vpstd, alpha=0.5)

Legend(f[1,2], [p, v, vp],["mean p angle", "mean v angle", "mean angle difference vp"])
display(f)
save(joinpath(base_folder,"mean_angular_dynamics.pdf"),f)


f = Figure(size=(2000,500));
ax=Axis(f[1,1], ylabel = "mean vp angle zoomed", xlabel="t", title="mean vp angle zoomed")

vp=plot!(ax, save_tax,vpm,colormap=:viridis)
band!(ax, save_tax, vpm-vpstd, vpm+vpstd, alpha=0.5)
ylims!(-0.5,0.5)

Legend(f[1,2], [vp],["mean angle difference vp"])
display(f)
save(joinpath(base_folder,"mean_vp_angle_zoomed.pdf"),f)

end


mean(vpm[end-10:end])


function construct_2d_rotation_matrix(angle)

    return [ cos(angle) -sin(angle); sin(angle) cos(angle) ]
end

Rot = construct_2d_rotation_matrix(-mean(vpm))

L = copy(D[:,:])
v0 = sa_file["frames"]["1"]["v0"][1]
J = sa_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]

for i in 1:Np

    L[2i-1:2i, 2i-1:2i]+=-J*v0/vrms[end] .*(I - Rot)

end


Leigenfact = eigen(L)





f = Figure(size=(4000,4000));
mm=1

for i in 1:4
    for j in 1:4
    print(mm)

    m = mm
    

    eigen_m_value = Leigenfact.values[m]
    eigen_m_vector =  reshape(Leigenfact.vectors[:,m], (2, round(Int64,length(Leigenfact.vectors[:,m])/2)) )
    eigen_m_x = eigen_m_vector[1,:]
    eigen_m_y = eigen_m_vector[2,:]


    ax=Axis(f[i,j], xlabel = "x", ylabel="y", aspect=1, title=string(eigen_m_value))
    scatter!(ax,x0,y0, markerspace=:data, markersize=2 .*R.^2, color = id,marker = Circle,alpha=0.1)

    arrows!(ax, x0, y0, real(eigen_m_x), real(eigen_m_y), lengthscale = 100)

    arrows!(ax, x0, y0, imag(eigen_m_x), imag(eigen_m_y), lengthscale = 100, color=:red)
    mm+=1
    end

end
display(f)


#Now do projections over time

Nframes = length(sa_file["frames"])
LNeigenvectors =  length(Leigenfact.values)
#Lpos_projections = zeros(Nframes,Neigenvectors )
Lvel_projections= zeros(Nframes,Neigenvectors )

x0interior =  extract_data_for_type("x", 1, reference_frame)
y0interior =  extract_data_for_type("y", 1, reference_frame)

vx0interior =  extract_data_for_type("vx", 1, reference_frame)
vy0interior =  extract_data_for_type("vy", 1, reference_frame)

velocity0interior = collect(Iterators.flatten(zip(vx0interior, vy0interior)))



@views for i in 1:Nframes

    #xi = extract_data_for_type("x", 1, sa_file["frames"][string(i)])
    #yi = extract_data_for_type("y", 1, sa_file["frames"][string(i)])

    vxi = extract_data_for_type("vx", 1, sa_file["frames"][string(i)])
    vyi = extract_data_for_type("vy", 1, sa_file["frames"][string(i)])

    #Interweave and calculate normalized displacement
    #displacement = collect(Iterators.flatten(zip(xi .- x0interior, yi .- y0interior)))

    #Interweave and calculate normalized velocity
    velocity = collect(Iterators.flatten(zip(vxi, vyi)))

    for j in 1:LNeigenvectors
        #Lpos_projections[i,j] = abs(sum( displacement .* Leigenfact.vectors[:,j] ))
        Lvel_projections[i,j] = real( (sum( velocity .* Leigenfact.vectors[:,j] ))/ (sum( velocity0interior .* Leigenfact.vectors[:,j])) )
    end

end

f = Figure();
ax=Axis(f[1,1], ylabel = "Re{<λL_m|v(t)>/<λL_m|v(0)>}", xlabel="t", title="L velocity projections")
save_tax = sa_file["integration_info"]["save_tax"]
Nplot = 100
#ylims!(1e-4,5e-1)
for i in  1:Nplot
    scatter!(ax, save_tax, Lvel_projections[:,i], color=abs(Leigenfact.values[i]), colorrange=(0,abs(Leigenfact.values[Nplot])), colormap=:viridis)
    lines!(ax, save_tax, exp.(-real(Leigenfact.values[i]).*save_tax).*cos.(imag(Leigenfact.values[i]).*save_tax), color=abs(Leigenfact.values[i]), colorrange=(0,abs(Leigenfact.values[Nplot])), colormap=:viridis)

end
display(f)
save(joinpath(figure_base_folder,"Lvelocity_projections.pdf"),f)


#animation of mode m 

for m=1:6
    eigen_m_value = Leigenfact.values[m]
    eigen_m_vector =  reshape(Leigenfact.vectors[:,m], (2, round(Int64,length(Leigenfact.vectors[:,m])/2)) )

    using GLMakie
    GLMakie.activate!(; focus_on_show=true)

    f = Figure(size=(500,500));
    t = Observable(0.)

    tax = 1:0.5:1000



    eigen_m_x = @lift( real(eigen_m_vector.* exp(- $t* eigen_m_value) )[1,:])
    eigen_m_y = @lift(real(eigen_m_vector.* exp(- $t* eigen_m_value) )[2,:])
    angles = @lift( angle.($eigen_m_x .+ $eigen_m_y .*1im) )

    ax=Axis(f[1,1], xlabel = "x", ylabel="y", aspect=1, title = @lift("t = $(round($t, digits = 1)), ev=  $(eigen_m_value)"))

    scatter!(ax,x0[Np+1:end],y0[Np+1:end], markerspace=:data, markersize=2 .*R[Np+1:end].^2, color = 0,marker = Circle,alpha=0.5,colormap=:hsv ,colorrange=(-pi,pi))
    scatter!(ax,x0[1:Np],y0[1:Np], markerspace=:data, markersize=2 .*R[1:Np].^2, color = angles,marker = Circle,alpha=0.5,colormap=:hsv ,colorrange=(-pi,pi))

    arrows!(ax, x0[1:Np], y0[1:Np], eigen_m_x, eigen_m_y, lengthscale = 100)
    display(f)

    record(f, "time_evolution_mode"*string(m)*".mp4", tax; visible=true, fps=120 ) do ti

        t[]=ti
    end
end


fig, ax,=plot(real(Leigenfact.values))
plot!(ax,imag(Leigenfact.values))



