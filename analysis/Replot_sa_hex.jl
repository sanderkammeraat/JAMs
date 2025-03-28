
include(joinpath("..","src","Engine.jl"))
using GLMakie

function main(base_folder)

    

    readdir(base_folder)

    Drfolders = readdir(base_folder)

    function make_movie(raw_data_file,save_folder)
        save_tax = raw_data_file["integration_info"]["save_tax"]

        J = raw_data_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]

        
        

        frame_numbers = 1:1:1000# length(save_tax)

        frames = raw_data_file["frames"]
        Dr = frames["1"]["Dr"][1]

        t = Observable(0.)


        x = Observable(frames["1"]["x"])
        y = Observable(frames["1"]["y"]) 

        vx = Observable(frames["1"]["vx"])
        vy = Observable(frames["1"]["vy"]) 

        px = Observable(frames["1"]["px"])
        py = Observable(frames["1"]["py"]) 
        type = frames["1"]["type"]

        t_stdp = Observable(Point2f[(0, std( px[][type.==1]) )])

        t_mp = Observable(Point2f[(0, mean( px[][type.==1]) )])

        R = Observable(frames["1"]["R"]) 
        
        id = frames["1"]["id"]
        Np = length(frames["1"]["x"])
        scaleup = maximum(frames["1"]["R"])+maximum(frames["1"]["x"])
        #Setup figure
        f = Figure(size=(1000,500));
        ax = Axis(f[1,1], aspect=DataAspect(),title = @lift("t = $(round($t, digits = 1)), Dr=$Dr, J=$J "), xlabel="x", ylabel="y");

        ax2 = Axis(f[1,2],title = @lift("t = $(round($t, digits = 1)), Dr=$Dr, J=$J "), xlabel="t", ylabel="std(px), mean(px)");


        #disks
        c = @lift( angle.($px+1im*$py) )
        s = @lift( 2. *$R.^2 )

        scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.1, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))



        #directors
        cp = @lift( angle.($px+1im*$py) )
        arrows!(ax, x,y, px,py , color=cp,  colormap=:hsv,colorrange=(-pi,pi))

        #velocity vectors
        cv = @lift( angle.($vx+1im*$vy) )
        arrows!(ax, x,y, vx,vy , color=cv,  colormap=:hsv,colorrange=(-pi,pi))

        plotpx = @lift(scaleup.*$(px)[type.==1]) 
        plotpy = @lift(scaleup.*$(py)[type.==1]) 

        scatter!(ax,plotpx,plotpy, color=cv,  colormap=:hsv,colorrange=(-pi,pi))


        lines!(ax2, t_stdp)

        lines!(ax2, t_mp)
        ylims!(ax2, (-1,1))



        display(f)

        record(f, joinpath(base_folder,save_folder,"simulation.mp4"), frame_numbers; visible=true) do i 

            stri = string(i)
            t[] = frames[stri]["t"]
            #tplot[]=push!(tplot[], frames[stri]["t"] )

            

            x[] = frames[stri]["x"]
            y[] = frames[stri]["y"]

            vx[] = frames[stri]["vx"]
            vy[] = frames[stri]["vy"]

            px[] = frames[stri]["px"]
            py[] = frames[stri]["py"]

            t_stdp[] = push!(t_stdp[], Point2f(frames[stri]["t"], std( px[][type.==1]) ) )

            t_mp[] = push!(t_mp[], Point2f(frames[stri]["t"], mean( px[][type.==1]) ) )

            xlims!(ax2, (0,1e-3+frames[stri]["t"]))
            
        end
    end



    @showprogress dt = 1 desc="Animating in progress..." showspeed=true for Drfolder in Drfolders

        folder_path = joinpath(base_folder, Drfolder)
        raw_data_file = jldopen( joinpath(folder_path, "raw_data.jld2"), "r" )

        make_movie(raw_data_file,folder_path)
        close(raw_data_file)

    end
end
base_folder = joinpath("/data1","kammeraat", "sa", "varyDr","J_1")
main(base_folder)

