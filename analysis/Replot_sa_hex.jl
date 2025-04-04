
include(joinpath("..","src","Engine.jl"))
include("AnalysisFunctions.jl")
include("AnalysisPipeline.jl")
using GLMakie
GLMakie.activate!()

function make_movie(raw_data_file,save_folder)
    save_tax = raw_data_file["integration_info"]["save_tax"]

    J = raw_data_file["system"]["forces"]["external_forces"]["self_align_with_v_unit_force"]["β"]

    frame_numbers = 1:1:length(save_tax)

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
    f = Figure(size=(2000,2000));
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

    vp_x = @lift( cos.($cp-$cv))

    vp_y = @lift( sin.($cp-$cv) )

    arrows!(ax, x,y, vx,vy , color=cv,  colormap=:hsv,colorrange=(-pi,pi))


    #arrows!(ax, x,y, vp_x,vp_y , color=cv,  colormap=:hsv,colorrange=(-pi,pi))

    plotpx = @lift(scaleup.*$(px)[type.==1]) 
    plotpy = @lift(scaleup.*$(py)[type.==1]) 

    scatter!(ax,plotpx,plotpy, color=cp,  colormap=:hsv,colorrange=(-pi,pi))


    lines!(ax2, t_stdp)

    lines!(ax2, t_mp)

    colsize!(f.layout, 1, Relative(3/4))
    ylims!(ax2, (-1,1))



    display(f)

    record(f, joinpath(save_folder,"Dr_$(Dr)_J_$(J).mp4"), frame_numbers; visible=true) do i 

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

function main(base_folder, animation_base_folder; raw_data_file_name="raw_data.jld2")

    tree = construct_folder_tree_param_param_seed(base_folder)

    for (param1, subdict) in tree
    
        @showprogress for (param2, seeddict) in subdict

            print(param1)
            #if param1=="Dr_1.0" && param2=="J_0.0"
    
                for (seed, seedpath) in  seeddict
        
                        raw_data_file_path = joinpath(seedpath, raw_data_file_name)
        
                        save_folder = joinpath( mkpath(joinpath(animation_base_folder,param1, param2,seed)))
        
                        raw_data_file = jldopen( raw_data_file_path, "r" )

                        make_movie(raw_data_file,save_folder)
                        close(raw_data_file)    
                end
            #end
        end
    end 
end

base_folder = joinpath("/data1","kammeraat", "sa", "phi_1","Nlin_20","vary_J_Dr")

#base_folder = joinpath(homedir(), "sa", "vary_J_Dr_largeN","simdata")
animation_base_folder = joinpath(base_folder,"movies")
main(joinpath(base_folder,"simdata"), animation_base_folder)

