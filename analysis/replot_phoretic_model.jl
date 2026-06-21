
include(joinpath("..","src","Engine.jl"))
include("AnalysisFunctions.jl")
include("AnalysisPipeline.jl")
using GLMakie
GLMakie.activate!()

@views function make_movie(raw_data_file,save_folder)
    save_tax = raw_data_file["integration_info"]["save_tax"]

    frame_numbers = 1:10:length(save_tax)

    um = 5
    scaleup=10

    frames = raw_data_file["frames"]

    Dr = frames["1"]["Dr"][1]

    t = Observable(0.)

    field_C = Observable(frames["1"]["field_C"])
    field_x_centers = frames["1"]["field_bin_centers_x"]
    field_y_centers = frames["1"]["field_bin_centers_y"]
    x = Observable(frames["1"]["x"].*um)
    y = Observable(frames["1"]["y"].*um) 

    xt = Observable(frames["1"]["x"].*um)
    yt = Observable(frames["1"]["y"].*um) 

    vx = Observable(frames["1"]["vx"].*um)
    vy = Observable(frames["1"]["vy"].*um) 

    px = Observable(frames["1"]["px"]*scaleup)
    py = Observable(frames["1"]["py"]*scaleup) 
    type = frames["1"]["type"]

    R =Observable(frames["1"]["R"]) 
    
    #Setup figure
    f = Figure(size=(1000,1000))#Figure(size=(2000,2000));
    ax = Axis(f[1,1], aspect=DataAspect(),title = @lift("t = $(round($t, digits = 1)) s"), xlabel="x (μm)", ylabel="y (μm)", xgridvisible=false, ygridvisible=false, xticks=[-750,-500,-250,0, 250, 500, 750], yticks=[-750,-500,-250,0, 250, 500, 750]);
    xlims!(ax, low=-750, high=750)
    ylims!(ax, low=-750, high=750)


    #disks
    c = type
    s = @lift( 2. *$R.*um)
    st = @lift( 2. *$R.*um/5)

    #heatmap!(ax,field_x_centers.*um, field_y_centers.*um, field_C, alpha=0.2,colormap=:viridis,colorrange=(0,1))

    cmap = :seismic
    heatmap!(ax,field_x_centers.*um,field_y_centers.*um,field_C, alpha=0.2,colormap=cmap,colorrange=(0.5,1.5))
    Colorbar(f[1,2], limits = (0.5, 1.5), label="Concentration c",colormap=(cmap,0.2))

    #vlines!(ax,field_x_centers, color="black", alpha=0.2)
    #hlines!(ax,field_y_centers, color="black", alpha=0.2)
    #scatter!(ax,xt,yt, color=:grey, alpha=0.1, markersize=2)
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1, colormap=Reverse(:seismic))

    



    #directors
    cp = @lift( angle.($px+1im*$py) )

    arrows2d!(ax, x,y, px,py , color=cp,  colormap=:hsv,colorrange=(-pi,pi))

    #velocity vectors
    cv = @lift( angle.($vx+1im*$vy) )

    vp_x = @lift( cos.($cp-$cv))

    vp_y = @lift( sin.($cp-$cv) )

    #arrows!(ax, x,y, vx,vy , color=cv,  colormap=:hsv,colorrange=(-pi,pi))
    delay=200
    display(f)
    save_path = mkpath( joinpath(save_folder,"simulation.mp4"))
    record(f,save_path, frame_numbers; visible=true) do i 

        stri = string(i)
        t[] = frames[stri]["t"]
        #tplot[]=push!(tplot[], frames[stri]["t"] )

        x[] = frames[stri]["x"].*um
        y[] = frames[stri]["y"].*um

 
        xt[] = vcat(xt[], frames[stri]["x"].*um)
        yt[] = vcat(yt[],frames[stri]["y"].*um)

        vx[] = frames[stri]["vx"].*um
        vy[] = frames[stri]["vy"].*um
        
        px[] = frames[stri]["px"]*scaleup
        py[] = frames[stri]["py"]*scaleup

        field_C[] = frames[stri]["field_C"]

        
    end
    close(raw_data_file)
end

base_folder = "/Users/kammeraat/Downloads/for_inference_v21_exp_rep/phi_0p01"

raw_data_file = jldopen(joinpath(base_folder,"simdata","raw_data.h5"),"r",iotype=IOStream )
make_movie(raw_data_file,joinpath(base_folder,"movies"))

close(raw_data_file)