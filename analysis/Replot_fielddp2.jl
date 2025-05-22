
include(joinpath("..","src","Engine.jl"))
include("AnalysisFunctions.jl")
include("AnalysisPipeline.jl")
using GLMakie
GLMakie.activate!()

@views function make_movie(raw_data_file,save_folder)
    save_tax = raw_data_file["integration_info"]["save_tax"]

    frame_numbers = 1:1:length(save_tax)

    frames = raw_data_file["frames"]

    Dr = frames["1"]["Dr"][1]

    t = Observable(0.)

    field_C = Observable(frames["1"]["field_C"][1])
    field_x_centers = frames["1"]["field_bin_centers"][1][1]
    field_y_centers = frames["1"]["field_bin_centers"][1][2]
    x = Observable(frames["1"]["x"])
    y = Observable(frames["1"]["y"]) 

    vx = Observable(frames["1"]["vx"])
    vy = Observable(frames["1"]["vy"]) 

    px = Observable(frames["1"]["px"])
    py = Observable(frames["1"]["py"]) 
    type = frames["1"]["type"]

    R = Observable(frames["1"]["R"]) 
    
    #Setup figure
    f = Figure(size=(1000,1000))#Figure(size=(2000,2000));
    ax = Axis(f[1,1], aspect=DataAspect(),title = @lift("t = $(round($t, digits = 1))"), xlabel="x", ylabel="y");
    xlims!(ax, low=field_x_centers[2], high=field_x_centers[end-1])
    ylims!(ax, low=field_y_centers[2], high=field_y_centers[end-1])


    #disks
    c = type
    s = @lift( 2. *$R)

    heatmap!(ax,field_x_centers, field_y_centers, field_C, alpha=0.2,colormap=:viridis,colorrange=(0,1))

    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1, colormap=Reverse(:seismic))



    #directors
    cp = @lift( angle.($px+1im*$py) )

    


    arrows!(ax, x,y, px,py , color=cp,  colormap=:hsv,colorrange=(-pi,pi))

    #velocity vectors
    cv = @lift( angle.($vx+1im*$vy) )

    vp_x = @lift( cos.($cp-$cv))

    vp_y = @lift( sin.($cp-$cv) )

    #arrows!(ax, x,y, vx,vy , color=cv,  colormap=:hsv,colorrange=(-pi,pi))



    display(f)
    save_path = mkpath( joinpath(save_folder,"anti-align_to_chemical_gradient.mp4"))
    record(f,save_path, frame_numbers; visible=true) do i 

        stri = string(i)
        t[] = frames[stri]["t"]
        #tplot[]=push!(tplot[], frames[stri]["t"] )

        

        x[] = frames[stri]["x"]
        y[] = frames[stri]["y"]

        vx[] = frames[stri]["vx"]
        vy[] = frames[stri]["vy"]

        px[] = frames[stri]["px"]
        py[] = frames[stri]["py"]

        field_C[] = frames[stri]["field_C"][1]

        
    end
    close(raw_data_file)
end

raw_data_file=jldopen(joinpath(homedir(),"ADCA","simulations","anti-align","simdata","raw_data.jld2"))




make_movie(raw_data_file,joinpath(homedir(),"ADCA","simulations","anti-align","movies"))
