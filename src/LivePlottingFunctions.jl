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
    heatmap!(ax,field_centers1,field_centers2,field_C, alpha=0.2,colormap=:viridis,colorrange=(0,1))
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



function plot_sized_points!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])

    

    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])

        R = @lift([p_i.R[1] for p_i in $cpsO])

        meshscatter!(ax,x,y,z, color=c, markersize =R, transparency=true)


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
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1, colormap=Reverse(:seismic))

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
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colorrange=(-0.05,0.05))

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



function setup_system_plotting(system_sizes,plot_functions,plotdim ,cpsO,cfsO,tO,fps;res=nothing)
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

        ax = Axis3(f[1, 1], xlabel = "x", ylabel="y", zlabel="z",  aspect = (1,system_sizes[2]/system_sizes[1],system_sizes[3]/system_sizes[1]), title=title)
        xlims!(ax,  -system_sizes[1]/2, system_sizes[1]/2)
        ylims!(ax, -system_sizes[2]/2, system_sizes[2]/2)
        zlims!(ax,  -system_sizes[3]/2, system_sizes[3]/2)

    end
    
    for plot_function in plot_functions
       ax=plot_function(f,ax,cpsO,cfsO)
    end
    display(f)
    return f, ax
end


