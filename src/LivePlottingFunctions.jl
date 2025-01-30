#Set default plotting backend
using GLMakie


function sphere!(x,y,z, R)

    φ = 0:π/100:2π

    θ = 0:π/200:π


    x.+= [R*cos(θ) * sin(φ) for θ in θ, φ in φ]

    y.+= [R*sin(θ)*sin(φ) for θ in θ, φ in φ]
    
    z.+= [R*cos(φ) for θ in θ, φ in φ]
    
    return x,y,z
end

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

function plot_points_on_plane!(ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])

    scatter!(ax,x,y, color=c)
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
        print(length(cpsO[][1].x))
        
        z = @lift([p_i.x[3] for p_i in $cpsO])

        R = @lift([p_i.a  for p_i in $cpsO])

        meshscatter!(ax,x,y,z, color=c, markersize =R, transparency=true)


    else
        s = @lift([2*p_i.a^2  for p_i in $cpsO])
        scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    end
    return ax

end


function new_plot_sized_points!(ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])

    

    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])

        R = @lift([p_i.R for p_i in $cpsO])

        meshscatter!(ax,x,y,z, color=c, markersize =R, transparency=true)


    else
        s = [2*p_i.R^2  for p_i in current_particle_state]
        scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    end

    return ax
end

function plot_type_sized_points!(ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.type for p_i in $cpsO])

    

    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])

        R = @lift([p_i.a for p_i in $cpsO])

        meshscatter!(ax,x,y,z, color=c, markersize =R, transparency=true)


    else
        s = @lift([2*p_i.a^2  for p_i in $cpsO])
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
        c = @lift([ p_i.p[1] for p_i in $cpsO])
        arrows!(ax, x, p , color=c,  colormap=:seismic,colorrange=(-1, 1))
        


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
        c = @lift([ norm(p_i.v) for p_i in $cpsO])
        arrows!(ax, x, v , color=c,  colormap=:seismic)
        


    else
        vx = @lift([ p_i.v[1] for p_i in $cpsO])
        vy = @lift([ p_i.v[2] for p_i in $cpsO])
        c = @lift([ norm(p_i.v) for p_i in $cpsO])
        arrows!(ax, x, y, vx, vy, color=c,  colormap=:plasma)


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