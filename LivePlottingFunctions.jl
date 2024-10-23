
function disk(x, y, R)
    
    ϕ = range(0, 2 * pi,100)
    return x .+ R*cos.(ϕ), y .+ R*sin.(ϕ)
end


#Define functions for live plotting
function plot_disks!(p, current_state)
        x = [p_i.x[1] for p_i in current_state]
        y = [p_i.x[2] for p_i in current_state]
        a = [ p_i.a  for p_i in current_state]
        plot!( disk.(x, y , a), seriestype=[:shape], lw=0.5, c=:blue, linecolor = :black, legend=false; fillalpha=0.2, aspect_ratio=1)
end

function plot_Swarmalators!(p, current_state)

    x = [p_i.x[1] for p_i in current_state]
    y = [p_i.x[2] for p_i in current_state]

    c = angle2range.([ p_i.θ[1] for p_i in current_state])

    scatter!(x,y , zcolor=c, colormap=:hsv, legend=false, colorbar = :none)

end

function plot_points!(p, current_state)

    x = [p_i.x[1] for p_i in current_state]
    y = [p_i.x[2] for p_i in current_state]

    c = [ p_i.id[1] for p_i in current_state]

    scatter!(x,y , zcolor=c, colormap=:hawaii, legend=false, colorbar = :none)

end

function plot_directors!(p, current_state)

    x = [p_i.x[1] for p_i in current_state]
    y = [p_i.x[2] for p_i in current_state]

    nx = cos.([ p_i.θ[1] for p_i in current_state])
    ny = sin.([ p_i.θ[1] for p_i in current_state])
    c = angle2range.([ p_i.θ[1] for p_i in current_state])
    c = [c c]'
    quiver!(x,y, quiver=(nx, ny),  line_z=repeat([c...], inner=2), colormap=:hsv,legend=false, clims=(0, 2*pi))

end

function plot_polar_disks!(p, current_state)

    x = [p_i.x[1] for p_i in current_state]
    y = [p_i.x[2] for p_i in current_state]
    a = [ p_i.a  for p_i in current_state]
    

    nx = cos.([ p_i.θ[1] for p_i in current_state])
    ny = sin.([ p_i.θ[1] for p_i in current_state])
    c = angle2range.([ p_i.θ[1] for p_i in current_state])
    c = [c c]'
    plot!( disk.(x, y , a), seriestype=[:shape], lw=0.5, c=:blue, linecolor = :black, legend=false; fillalpha=0.2, aspect_ratio=1)

    quiver!(x,y, quiver=(nx, ny),  line_z=repeat([c...], inner=2), colormap=:hsv,legend=false, clims=(0, 2*pi))

end

function angle2range(angle)
    if angle>=0
        return angle % (2*pi)
    else
        return angle % (2*pi) + (2*pi)
    end

end