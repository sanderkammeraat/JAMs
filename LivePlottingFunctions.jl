
function disk(x, y, R)
    
    ϕ = range(0, 2 * pi,100)
    return x .+ R*cos.(ϕ), y .+ R*sin.(ϕ)
      
end


#Define functions for live plotting
function plot_disks!(p, current_state, n, Tplot)
    if n%Tplot==0
        
        empty!(p)

        for p_i in current_state
            plot!( disk(p_i.x[1], p_i.x[2] , p_i.a), seriestype=[:shape], lw=0.5, c=:blue, linecolor = :black, legend=false; fillalpha=0.2, aspect_ratio=1)
        end

        display(p)
    end
end

function plot_Swarmalators!(p, current_state, n, Tplot)


    if n%Tplot==0
        
        empty!(p)

        x = [p_i.x[1] for p_i in current_state]
        y = [p_i.x[2] for p_i in current_state]

        c = [ (p_i.θ[1] % (2*pi) + 2 *pi) for p_i in current_state]

        scatter!(x,y , zcolor=c, colormap=:hsv, legend=false, colorbar = :none)

        display(p)
    end
end

function plot_points!(p, current_state, n, Tplot)


    if n%Tplot==0
        
        empty!(p)

        x = [p_i.x[1] for p_i in current_state]
        y = [p_i.x[2] for p_i in current_state]

        c = [ p_i.id[1] for p_i in current_state]

        scatter!(x,y , zcolor=c, colormap=:hawaii, legend=false, colorbar = :none)

        display(p)
    end
end