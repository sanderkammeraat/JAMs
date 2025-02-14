#Set default plotting backend
using GLMakie


function plot_points!(ax, cpsO, cfsO)


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

function plot_field_magnitude!(ax, cpsO, cfsO)
    
    field_centers1 = @lift($(cfsO)[1].bin_centers[1])
    field_centers2 = @lift($(cfsO)[1].bin_centers[2])
    field_C = @lift($(cfsO)[1].C)
    heatmap!(ax,field_centers1,field_centers2,field_C, alpha=0.2,colormap=:viridis,colorrange=(0,0.4))
    return ax
end


function plot_type_points!(ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.type for p_i in $cpsO])
    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])
        meshscatter!(ax,x,y,z, color=c)

    else
        scatter!(ax,x,y, color=c)

    end
    return ax
end


function plot_Swarmalators!(ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift(angle2range.([ p_i.ϕ[1] for p_i in $cpsO]))

    scatter!(ax,x,y , color=c, colormap=:hsv, colorrange=(0, 2*pi))
    return ax
end



function plot_sized_points!(ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])

    

    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])

        R = @lift([p_i.R for p_i in $cpsO])

        meshscatter!(ax,x,y,z, color=c, markersize =R, transparency=true)


    else
        s = @lift([2*p_i.R^2  for p_i in $cpsO])
        scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    end

    return ax
end



function plot_disks!(ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R^2  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    return ax
end

function plot_disks_vx!(ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.v[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R^2  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colorrange=(-0.05,0.05))

    return ax
end


function plot_type_sized_points!(ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.type for p_i in $cpsO])

    

    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])

        R = @lift([p_i.R for p_i in $cpsO])

        meshscatter!(ax,x,y,z, color=c, markersize =R, transparency=true)


    else
        s = @lift([2*p_i.R^2  for p_i in $cpsO])
        scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    end
    return ax

end

function plot_directors!(ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    
    if length(cpsO[][1].x)>2

        x = @lift([ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in $cpsO])
        p = @lift([ Point3f( p_i.p[1],p_i.p[2],p_i.p[3]) for p_i in $cpsO])
        c = @lift([ angle(p_i.p[1]+1im*p_i.p[2]) for p_i in $cpsO])
        arrows!(ax, x, p , color=c,  colormap=:hsv,colorrange=(-pi,pi))
        


    else
        nx = @lift(cos.([ p_i.θ[1] for p_i in $cpsO]))
        ny = @lift(sin.([ p_i.θ[1] for p_i in $cpsO]))
        c = @lift(angle2range.([ p_i.θ[1] for p_i in $cpsO]))
        arrows!(ax, x, y, nx, ny, color=c,  colormap=:hsv,colorrange=(0, 2*pi))


    end
    return ax

end

function plot_velocity_vectors!(ax,cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    
    if length(cpsO[][1].x)>2

        x = @lift([ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in $cpsO])
        v = @lift([ Point3f( p_i.v[1],p_i.v[2],p_i.v[3]) for p_i in $cpsO])
        c = @lift([ angle(p_i.v[1]+1im*p_i.v[2]) for p_i in $cpsO])
        arrows!(ax, x, v , color=c,  colormap=:hsv,colorrange=(-pi,pi))
        


    else
        vx = @lift([ p_i.v[1] for p_i in $cpsO])
        vy = @lift([ p_i.v[2] for p_i in $cpsO])
        c = @lift([ angle(p_i.v[1]+1im*p_i.v[2]) for p_i in $cpsO])
        arrows!(ax, x, y, vx, vy, color=c, colormap=:hsv,colorrange=(-pi,pi))


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



function setup_system_plotting(system_sizes,plot_functions,plotdim ,cpsO,cfsO,tO,res=nothing)
    GLMakie.activate!()
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
       ax=plot_function(ax,cpsO,cfsO)
    end
    display(f)
    return f, ax
end



function make_movie(SIM, folder_path, file_name, plot_functions, fps,plotdim=nothing)

    mkpath(save_path)

    cpsO = Observable(SIM.particle_states[1])

    cfsO = Observable(SIM.field_states[1])

    tO = Observable(SIM.tsax[1])

    t_indices = range(1,length(SIM.tsax))

    fig, ax = setup_system_plotting(SIM.system.sizes,plot_functions,plotdim ,cpsO,cfsO,tO,(500,500))
    
    record(fig, folder_path*file_name, t_indices; framerate=fps, compression=30) do t_index 
        cpsO[] = SIM.particle_states[t_index]
        cfsO[] = SIM.field_states[t_index]
        tO[] = SIM.tsax[t_index]
    end


end

function make_snapshot(SIM, folder_path, file_name, plot_functions, frame_index,plotdim=nothing)

    mkpath(folder_path)

    cpsO = Observable(SIM.particle_states[frame_index])

    cfsO = Observable(SIM.field_states[frame_index])

    tO = Observable(SIM.tsax[frame_index])

    t_indices = range(1,length(SIM.tsax))

    fig, ax = setup_system_plotting(SIM.system.sizes,plot_functions,plotdim ,cpsO,cfsO,tO,(1000,1000))

    file_path = folder_path*file_name

    save(file_path,fig)



end
