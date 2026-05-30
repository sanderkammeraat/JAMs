#Set default plotting backend
using GLMakie

if pkgversion(Makie) >= v"0.23-"
    gen_arrows_2d!(args...; kwargs...) = arrows2d!(args...; kwargs...)
    gen_arrows_3d!(args...; kwargs...) = arrows3d!(args...; kwargs...)
else
    gen_arrows_2d!(args...; kwargs...) = arrows!(args...; kwargs...)
    gen_arrows_3d!(args...; kwargs...) = arrows!(args...; kwargs...)
end

function plot_points!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])
    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])
        meshscatter!(ax,x,y,z, color=c)

    else
        scatter!(ax,x,y, color=c)

    end
    return ax
end

#Experimental
function plot_trajectories!(f,ax, cpsO, cfsO)
    maxlen = 100
    N = length(cpsO[])

    trajectories = fill(Point3f(NaN, NaN, NaN), (maxlen+1) * N)
    indices = zeros(Int, N)
    
    obs = Observable(trajectories)
    color_idx =repeat(1:N, inner=maxlen+1)
    lines!(ax, obs, color=color_idx, colorrange=(1, N), colormap=:viridis)
    
    on(cpsO) do particles
        for i in 1:N
            p = particles[i]
            idx = (indices[i] % maxlen) + 1
            indices[i] = idx
            
            base_idx = (i - 1) * (maxlen+1)
            trajectories[base_idx + idx] = Point3f(p.x[1], p.x[2], p.x[3])
            if idx <maxlen
                trajectories[base_idx + idx+ 1] = Point3f(NaN, NaN, NaN)
            else
                trajectories[base_idx + 2] = Point3f(NaN, NaN, NaN)
            end
    
            
        end
        notify(obs)
    end
    
    return ax
end

function plot_director_points!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])
    

    c = @lift([ p_i.id[1] for p_i in $cpsO])
    if length(cpsO[][1].x)>2


        
        z = @lift([p_i.x[3] for p_i in $cpsO])
        xpx = @lift([p_i.p[1]+p_i.x[1] for p_i in $cpsO])
        ypy = @lift([p_i.p[2]+p_i.x[2] for p_i in $cpsO])

        meshscatter!(ax,xpx,ypy,z, color=:red)

    else
        scatter!(ax,x,y, color=c)

    end
    return ax
end

function plot_field_magnitude!(f,ax, cpsO, cfsO)
    
    field_centers1 = @lift($(cfsO)[1].bin_centers[1])
    field_centers2 = @lift($(cfsO)[1].bin_centers[2])
    
    field_C = @lift($(cfsO)[1].C)
    heatmap!(ax,field_centers1,field_centers2,field_C, alpha=0.2,colormap=:jet,colorrange=(0.0,1))
    Colorbar(f[1,2], limits = (0.0, 1), label="Concentration c",colormap=(:jet,0.2))
    return ax
end

function plot_field_log_magnitude!(f,ax, cpsO, cfsO)
    
    field_centers1 = @lift($(cfsO)[1].bin_centers[1])
    field_centers2 = @lift($(cfsO)[1].bin_centers[2])
    
    field_C = @lift(log10.($(cfsO)[1].C))
    heatmap!(ax,field_centers1,field_centers2,field_C, alpha=0.2,colormap=:jet,colorrange=(-1,0))
    Colorbar(f[1,2], limits = (-1, 0), label="Concentration c",colormap=:jet)
    return ax
end



function plot_potential!(f,ax, cpsO, cfsO)
    
    field_centers1 = @lift($(cfsO)[1].bin_centers[1])
    field_centers2 = @lift($(cfsO)[1].bin_centers[2])
    field_C = @lift($(cfsO)[1].C)
    Cmax = maximum(cfsO[][1].C)
    Cmin = minimum(cfsO[][1].C)
    surface!(ax,field_centers1,field_centers2,field_C,colormap=Reverse(:gist_rainbow), alpha=0.5, colorrange=(Cmin, Cmax), color=field_C,shading=NoShading)
    Colorbar(f[1,2], limits = (Cmin, Cmax), label="Pot. energy U",colormap=Reverse(:gist_rainbow))
    return ax
end

function plot_field_magnitude_wgrid!(f,ax, cpsO, cfsO)
    
    field_centers1 = @lift($(cfsO)[1].bin_centers[1])
    field_centers2 = @lift($(cfsO)[1].bin_centers[2])
    field_C = @lift($(cfsO)[1].C)
    heatmap!(ax,field_centers1,field_centers2,field_C, alpha=0.2,colormap=:viridis,colorrange=(0,1))

    vlines!(ax,field_centers1, color="black", alpha=0.2)
    hlines!(ax,field_centers2, color="black", alpha=0.2)
    return ax
end

function plot_type_points!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.type[1] for p_i in $cpsO])
    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])
        meshscatter!(ax,x,y,z, color=c)

    else
        scatter!(ax,x,y, color=c)

    end
    return ax
end

function plot_type_points!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.type[1] for p_i in $cpsO])
    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])
        meshscatter!(ax,x,y,z, color=c)

    else
        scatter!(ax,x,y, color=c)

    end
    return ax
end

function plot_shape_disks!(f,ax, cpsO, cfsO)
    c = @lift([ p_i.id[1] for p_i in $cpsO])
    for j = 1:size(cpsO[][1].xe)[1]

        xes = @lift([p_i.xe[j,1] for p_i in $cpsO] )

        yes = @lift([p_i.xe[j,2] for p_i in $cpsO] )

        s = @lift([2*p_i.re[j] for p_i in $cpsO])
    
        scatter!(ax,xes,yes, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.5, strokecolor=:black, strokewidth=1)
    end
    return ax
end

function plot_shape_disks_orientation!(f,ax, cpsO, cfsO)
    c = @lift([ angle(p_i.p[1]+1im*p_i.p[2]) for p_i in $cpsO])
    for j = 1:size(cpsO[][1].xe)[1]

        xes = @lift([p_i.xe[j,1] for p_i in $cpsO] )

        yes = @lift([p_i.xe[j,2] for p_i in $cpsO] )

        s = @lift([2*p_i.re[j] for p_i in $cpsO])
    
        scatter!(ax,xes,yes, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.5, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))
    end
    return ax
end

function plot_shape_disks_type!(f,ax, cpsO, cfsO)
    c = @lift([ p_i.type[1] for p_i in $cpsO])
    for j = 1:size(cpsO[][1].xe)[1]

        xes = @lift([p_i.xe[j,1] for p_i in $cpsO] )

        yes = @lift([p_i.xe[j,2] for p_i in $cpsO] )

        s = @lift([2*p_i.re[j] for p_i in $cpsO])
    
        scatter!(ax,xes,yes, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.5, strokecolor=:black, strokewidth=1,colormap=Reverse(:seismic))
    end
    return ax
end

function plot_shape_points!(f,ax, cpsO, cfsO)
    c = @lift([ p_i.id[1] for p_i in $cpsO])
    for j = 1:size(cpsO[][1].xe)[1]

        xes = @lift([p_i.xe[j,1] for p_i in $cpsO] )

        yes = @lift([p_i.xe[j,2] for p_i in $cpsO] )

        s = @lift([p_i.re[1] for p_i in $cpsO])
    
        scatter!(ax,xes,yes, color=c)
    end
    return ax
end




function plot_Swarmalators!(f,ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift(angle2range.([ p_i.ϕ[1] for p_i in $cpsO]))

    scatter!(ax,x,y , color=c, colormap=:hsv, colorrange=(0, 2*pi))
    return ax
end


function plot_sphere!(f,ax, cpsO, cfsO)


    R = norm(cpsO[][1].x)
    meshscatter!(ax,0,0,0, markersize=R, color=:white)

    return ax
end


function plot_sized_points!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])

    

    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])

        R = @lift([p_i.R[1] for p_i in $cpsO])

        meshscatter!(ax,x,y,z, color=c, markersize =R, transparency=true, colormap=(:viridis,0.2))


    else
        s = @lift([2*p_i.R[1]  for p_i in $cpsO])
        scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    end

    return ax
end


function plot_polymers!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.pol_id[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1, colormap=:rainbow1)

    return ax
end



function plot_disks!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    return ax
end



function ellipse(p_i)

    lmda_major = 0.5*(p_i.Lambda[1,1]+p_i.Lambda[2,2]) + sqrt(0.25*(p_i.Lambda[1,1]-p_i.Lambda[2,2])^2 + p_i.Lambda[1,2]^2)

    lmda_minor =0.5*(p_i.Lambda[1,1]+p_i.Lambda[2,2]) -sqrt(0.25*(p_i.Lambda[1,1]-p_i.Lambda[2,2])^2 + p_i.Lambda[1,2]^2)

    s = range(0, 2pi, length=40)

    x_e = lmda_major .* cos.(s)

    y_e = lmda_minor .* sin.(s)

   x_e_rot = p_i.p[1] * x_e  - p_i.p[2] * y_e
   y_e_rot= p_i.p[2] * x_e  + p_i.p[1] *  y_e

   return p_i.x[1] .+ x_e_rot, p_i.x[2] .+ y_e_rot

end



function plot_ellipses!(f,ax, cpsO, cfsO)

    xs = @lift([ ellipse( p_i )[1] for p_i in $cpsO])
    ys = @lift([ ellipse( p_i )[2] for p_i in $cpsO])


    for i=1:length(cpsO[])

        xsi = @lift( $xs[i])
        ysi = @lift($ys[i])

        poly!(ax, xsi, ysi, color=i, colorrange=(1,length(cpsO[])),alpha=0.7, strokecolor=:black, strokewidth=2)
    end

    return ax
end

function plot_transparant_disks!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])
    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=:white, markersize =s,marker = Circle, markerspace=:data,alpha=0.1, strokecolor=:black, strokewidth=.2)

    return ax
end



function plot_disks_type!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.type[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1, colormap=:RdYlGn)

    return ax
end


function plot_disks_uw!(f,ax, cpsO, cfsO)


    x = @lift([p_i.xuw[1] for p_i in $cpsO])
    y = @lift([p_i.xuw[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    return ax
end

function plot_disks_vx!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.v[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colorrange=(-0.01,0.01),colormap=:seismic)

    return ax
end

function plot_disks_orientation!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ angle(p_i.p[1]+1im*p_i.p[2]) for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))

    return ax
end

function plot_disks_v_orientation!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ angle(p_i.v[1]+1im*p_i.v[2]) for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))

    return ax
end

function plot_disks_nematic_orientation!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ angle( exp(angle(p_i.p[1]+1im*p_i.p[2])*2  * 1im)) for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))

    return ax
end

function plot_disks_vp_phase_difference!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ angle( exp(1im* (angle(p_i.v[1]+1im*p_i.v[2]) - angle(p_i.p[1]+1im*p_i.p[2])) ) )  for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))

    return ax
end


function plot_type_sized_points!(f,ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.type[1] for p_i in $cpsO])

    

    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])

        R = @lift([p_i.R for p_i in $cpsO])

        meshscatter!(ax,x,y,z, color=c, markersize =R, transparency=true,markerspace=:data)


    else
        s = @lift([2*p_i.R[1]  for p_i in $cpsO])
        scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    end
    return ax

end

function plot_directors!(f,ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    
    if length(cpsO[][1].x)>2

        x = @lift([ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in $cpsO])
        p = @lift([ Point3f( p_i.p[1],p_i.p[2],p_i.p[3]) for p_i in $cpsO])
        c = @lift([ angle(p_i.p[1]+1im*p_i.p[2]) for p_i in $cpsO])
        gen_arrows_3d!(ax, x, p , color=c,  colormap=:hsv,colorrange=(-pi,pi))
        


    else
        nx = @lift(cos.([ p_i.θ[1] for p_i in $cpsO]))
        ny = @lift(sin.([ p_i.θ[1] for p_i in $cpsO]))
        c = @lift(angle2range.([ p_i.θ[1] for p_i in $cpsO]))
        gen_arrows_2d!(ax, x, y, nx, ny, color=c,  colormap=:hsv,colorrange=(0, 2*pi))


    end
    return ax

end

function plot_nematic_directors!(f,ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    
    if length(cpsO[][1].x)>2

        x = @lift([ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in $cpsO])
        p = @lift([ Point3f( p_i.p[1],p_i.p[2],p_i.p[3]) for p_i in $cpsO])

        pmin = @lift(-1*[ Point3f( p_i.p[1],p_i.p[2],p_i.p[3]) for p_i in $cpsO])


        c = @lift([ angle( exp(angle(p_i.p[1]+1im*p_i.p[2])*2  * 1im)) for p_i in $cpsO])
        gen_arrows_3d!(ax, x, p , color=c,  colormap=:hsv,colorrange=(-pi,pi))
        gen_arrows_3d!(ax, x, pmin , color=c,  colormap=:hsv,colorrange=(-pi,pi))


    else
        nx = @lift(cos.([ p_i.θ[1] for p_i in $cpsO]))
        ny = @lift(sin.([ p_i.θ[1] for p_i in $cpsO]))


        nxmin = @lift(-1*cos.([ p_i.θ[1] for p_i in $cpsO]))
        nymin = @lift(-1*sin.([ p_i.θ[1] for p_i in $cpsO]))


        c = @lift([ angle( exp(angle(p_i.p[1]+1im*p_i.p[2])*2  * 1im)) for p_i in $cpsO])
        gen_arrows_2d!(ax, x, y, nx, ny, color=c,  colormap=:hsv,colorrange=(0, 2*pi))
        gen_arrows_2d!(ax, x, y, nxmin, nymin, color=c,  colormap=:hsv,colorrange=(0, 2*pi))

    end
    return ax

end

function plot_velocity_vectors!(f,ax,cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    
    if length(cpsO[][1].x)>2

        x = @lift([ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in $cpsO])
        v = @lift([ Point3f( p_i.v[1],p_i.v[2],p_i.v[3]) for p_i in $cpsO])
        c = @lift([ angle(p_i.v[1]+1im*p_i.v[2]) for p_i in $cpsO])
        gen_arrows_3d!(ax, x, v , color=c,  colormap=:hsv,colorrange=(-pi,pi))
        


    else
        vx = @lift([ p_i.v[1] for p_i in $cpsO])
        vy = @lift([ p_i.v[2] for p_i in $cpsO])
        c = @lift([ angle(p_i.v[1]+1im*p_i.v[2]) for p_i in $cpsO])
        gen_arrows_2d!(ax, x, y, vx, vy, color=c, colormap=:hsv,colorrange=(-pi,pi))


    end
    return ax
end



function angle2range(angle)
    if angle>=0
        return angle % (2*pi)
    else
        return angle % (2*pi) + (2*pi)
    end

end



# function setup_system_plotting(system_sizes,plot_functions,plotdim ,cpsO,cfsO,tO,fps;res=nothing)
#     GLMakie.activate!(; focus_on_show=true, title= "GLMakie: JAMs simulation", framerate=fps)
#     if !isnothing(res)
#         f = Figure(size=res)
#     else
#         f=Figure()
#     end
#     title = @lift("t = $($tO)")

#     if !isnothing(plotdim)
#         plotdim_set = plotdim
#     else
#         plotdim_set = length(system_sizes)

#     end

#     if plotdim_set==2
#         ax = Axis(f[1, 1], xlabel = "x", ylabel="y",  aspect =system_sizes[1]/system_sizes[2], title=title )
#         xlims!(ax, -system_sizes[1]/2, system_sizes[1]/2)
#         ylims!(ax,  -system_sizes[2]/2, system_sizes[2]/2)

#     elseif plotdim_set==3

#         ax = Axis3(f[1, 1], xlabel = "x", ylabel="y", zlabel="z",  aspect = (1,system_sizes[2]/system_sizes[1],system_sizes[3]/system_sizes[1]), title=title)
#         xlims!(ax,  -system_sizes[1]/2, system_sizes[1]/2)
#         ylims!(ax, -system_sizes[2]/2, system_sizes[2]/2)
#         zlims!(ax,  -system_sizes[3]/2, system_sizes[3]/2)

#     end
    
#     for plot_function in plot_functions
#        ax=plot_function(f,ax,cpsO,cfsO)
#     end
#     display(f)
#     return f, ax
# end

function setup_system_plotting(system_sizes,plot_functions,plotdim ,cpsO,cfsO,tO,fps;res=nothing,sbs=false)
    GLMakie.activate!(; focus_on_show=true, title= "GLMakie: JAMs simulation", framerate=fps)
    if !isnothing(res)
        f = Figure(size=res)
    else
        f=Figure()
    end
    title = @lift("t = $($tO)")

    if !isnothing(plotdim)
        plotdim_set = plotdim
    else
        plotdim_set = length(system_sizes)

    end

    if plotdim_set==2
        ax = Axis(f[1, 1], xlabel = "x", ylabel="y",  aspect =system_sizes[1]/system_sizes[2], title=title )
        xlims!(ax, -system_sizes[1]/2, system_sizes[1]/2)
        ylims!(ax,  -system_sizes[2]/2, system_sizes[2]/2)

    elseif plotdim_set==3
        if !sbs
        ax = Axis3(f[1, 1], xlabel = "x", ylabel="y", zlabel="z",  aspect = (1,system_sizes[2]/system_sizes[1],system_sizes[3]/system_sizes[1]), title=title)
        xlims!(ax,  -system_sizes[1]/2, system_sizes[1]/2)
        ylims!(ax, -system_sizes[2]/2, system_sizes[2]/2)
        zlims!(ax,  -system_sizes[3]/2, system_sizes[3]/2)

    else
        eye_separation=0.03f0
        ax_left  = Axis3(f[1, 1], xlabel="x", ylabel="y", zlabel="z",aspect=(1,system_sizes[2]/system_sizes[1],system_sizes[3]/system_sizes[1]), title=title)
        ax_right = Axis3(f[1, 2], xlabel="x", ylabel="y", zlabel="z", aspect=(1,system_sizes[2]/system_sizes[1],system_sizes[3]/system_sizes[1]), title=title)
        for ax in (ax_left, ax_right)
                    xlims!(ax, -system_sizes[1]/2, system_sizes[1]/2)
                    ylims!(ax, -system_sizes[2]/2, system_sizes[2]/2)
                    zlims!(ax, -system_sizes[3]/2, system_sizes[3]/2)
        end

        syncing = Ref(false)

        on(ax_left.finallimits) do new_limits
            syncing[] && return
            syncing[] = true
            # Update right axis with the new limits
            xlims!(ax_right, new_limits.origin[1], new_limits.origin[1] + new_limits.widths[1])
            ylims!(ax_right, new_limits.origin[2], new_limits.origin[2] + new_limits.widths[2])
            zlims!(ax_right, new_limits.origin[3], new_limits.origin[3] + new_limits.widths[3])
            syncing[] = false
        end
        
        on(ax_right.finallimits) do new_limits
            syncing[] && return
            syncing[] = true
            # Update left axis with the new limits
            xlims!(ax_left, new_limits.origin[1], new_limits.origin[1] + new_limits.widths[1])
            ylims!(ax_left, new_limits.origin[2], new_limits.origin[2] + new_limits.widths[2])
            zlims!(ax_left, new_limits.origin[3], new_limits.origin[3] + new_limits.widths[3])
            syncing[] = false
        end

        on(ax_left.azimuth) do az
            syncing[] && return
            syncing[] = true
            ax_right.azimuth[] = az + 2f0 * eye_separation  # Reversed sign
            syncing[] = false
        end
        on(ax_right.azimuth) do az
            syncing[] && return
            syncing[] = true
            ax_left.azimuth[] = az - 2f0 * eye_separation  # Reversed sign
            syncing[] = false
        end
        # ── Elevation (identical for both eyes) ───────────────────────────────
        on(ax_left.elevation) do el
            syncing[] && return
            syncing[] = true
            ax_right.elevation[] = el
            syncing[] = false
        end
        on(ax_right.elevation) do el
            syncing[] && return
            syncing[] = true
            ax_left.elevation[] = el
            syncing[] = false
        end
        ax_right.azimuth[] = ax_left.azimuth[] + 2f0 * eye_separation
        end

    end
    if !sbs || plotdim_set==2
        for plot_function in plot_functions
            ax=plot_function(f,ax,cpsO,cfsO)
        end
        display(f)
        return f, ax
    else
        for plot_function in plot_functions
            ax_left=plot_function(f,ax_left,cpsO,cfsO)
            ax_right=plot_function(f,ax_right,cpsO,cfsO)
        end
        display(f)
        return f, (ax_left, ax_right)
    end

end




