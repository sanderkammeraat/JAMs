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

function plot_points!(ax, current_state)


    x = [p_i.x[1] for p_i in current_state]
    y = [p_i.x[2] for p_i in current_state]
    c = [ p_i.id[1] for p_i in current_state]
    if length(current_state[1].x)>2
        
        z = [p_i.x[3] for p_i in current_state]
        meshscatter!(ax,x,y,z, color=c)

    else
        scatter!(ax,x,y, color=c)

    end

end


function plot_Swarmalators!(ax, current_state)

    x = [p_i.x[1] for p_i in current_state]
    y = [p_i.x[2] for p_i in current_state]

    c = angle2range.([ p_i.ϕ[1] for p_i in current_state])

    scatter!(ax,x,y , color=c, colormap=:hsv, colorrange=(0, 2*pi))

end



function plot_sized_points!(ax, current_state)

    x = [p_i.x[1] for p_i in current_state]
    y = [p_i.x[2] for p_i in current_state]

    c = [ p_i.id[1] for p_i in current_state]

    

    if length(current_state[1].x)>2
        
        z = [p_i.x[3] for p_i in current_state]

        R = [p_i.a  for p_i in current_state]

        meshscatter!(ax,x,y,z, color=c, markersize =R, transparency=true)


    else
        s = [2*p_i.a^2  for p_i in current_state]
        scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    end


end


function plot_directors!(ax, current_state)

    x = [p_i.x[1] for p_i in current_state]
    y = [p_i.x[2] for p_i in current_state]

    
    if length(current_state[1].x)>2

        x = [ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in current_state]
        p = [ Point3f( p_i.p[1],p_i.p[2],p_i.p[3]) for p_i in current_state]
        c = [ p_i.p[1] for p_i in current_state]
        arrows!(ax, x, p , color=c,  colormap=:seismic,colorrange=(-1, 1))
        


    else
        nx = cos.([ p_i.θ[1] for p_i in current_state])
        ny = sin.([ p_i.θ[1] for p_i in current_state])
        c = angle2range.([ p_i.θ[1] for p_i in current_state])
        arrows!(ax, x, y, nx, ny, color=c,  colormap=:hsv,colorrange=(0, 2*pi))


    end


end

function plot_velocity_vectors!(ax, current_state)

    x = [p_i.x[1] for p_i in current_state]
    y = [p_i.x[2] for p_i in current_state]

    
    if length(current_state[1].x)>2

        x = [ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in current_state]
        v = [ Point3f( p_i.v[1],p_i.v[2],p_i.v[3]) for p_i in current_state]
        c = [ norm(p_i.v) for p_i in current_state]
        arrows!(ax, x, v , color=c,  colormap=:seismic)
        


    else
        vx = [ p_i.v[1] for p_i in current_state]
        vy = [ p_i.v[2] for p_i in current_state]
        c = [ norm(p_i.v) for p_i in current_state]
        arrows!(ax, x, y, vx, vy, color=c,  colormap=:plasma)


    end


end



function angle2range(angle)
    if angle>=0
        return angle % (2*pi)
    else
        return angle % (2*pi) + (2*pi)
    end

end