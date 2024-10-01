
#Define functions for live plotting
function plot_disks(current_state, n, Tplot)
    if n%Tplot==0
        
        x = [p_i.x[1] for p_i in current_state]
        y = [p_i.x[2] for p_i in current_state]
        s=[sqrt(2)*4*p_i.a/sqrt(system.sizes[1]*system.sizes[2]) for p_i in current_state]
        c = [p_i.id for p_i in current_state]

        
        #S=scatter(x,y, xlimits = (0,system.sizes[1]), ylimits=(0, system.sizes[2]), legend=false, zcolor=c, color=:hawaii, markersize=s,markerspace=SceneSpace ,aspect_ratio = :equal)
        f = Figure()

        ax = Axis(f[1, 1], xlabel = "x", ylabel = "y",aspect=DataAspect())
        
        scatter!(ax, x,y, markersize=s, color=c,markerspace=SceneSpace)

        xlims!(ax,(0, system.sizes[1]))
        ylims!(ax,(0, system.sizes[2]))
        display(f)
    end
end

function plot_Swarmalators(current_state, n, Tplot)


    if n%Tplot==0
        
        x = [p_i.x[1] for p_i in current_state]
        y = [p_i.x[2] for p_i in current_state]
        s=1 * sqrt(2) * 4/sqrt(system.sizes[1]*system.sizes[2])
        c = [ (p_i.θ[1] % (2*pi) + 2 *pi) for p_i in current_state]

        
        #S=scatter(x,y, xlimits = (0,system.sizes[1]), ylimits=(0, system.sizes[2]), legend=false, zcolor=c, color=:hawaii, markersize=s,markerspace=SceneSpace ,aspect_ratio = :equal)
        f = Figure()

        ax = Axis(f[1, 1], xlabel = "x", ylabel = "y",aspect=DataAspect())
        
        scatter!(ax, x,y, markersize=s, color=c, colormap=:hsv, markerspace=SceneSpace)

        xlims!(ax,(0, system.sizes[1]))
        ylims!(ax,(0, system.sizes[2]))
        display(f)
    end
end