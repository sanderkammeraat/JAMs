
include(joinpath("..","src","Engine.jl"))
include("AnalysisFunctions.jl")
include("AnalysisPipeline.jl")
using GLMakie
GLMakie.activate!()

@views function make_movie(raw_data_file,save_folder)
    save_tax = raw_data_file["integration_info"]["save_tax"]

    frame_numbers = 1:10:length(save_tax)

    um = 5

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

    px = Observable(frames["1"]["px"])
    py = Observable(frames["1"]["py"]) 
    type = frames["1"]["type"]

    R = Observable(frames["1"]["R"]) 
    
    #Setup figure
    f = Figure(size=(1000,1000))#Figure(size=(2000,2000));
    ax = Axis(f[1,1], aspect=DataAspect(),title = @lift("t = $(round($t, digits = 1))"), xlabel="x", ylabel="y");
    xlims!(ax, low=1*field_x_centers[2].*um, high=1*field_x_centers[end-1].*um)
    ylims!(ax, low=1*field_y_centers[2].*um, high=1*field_y_centers[end-1].*um)


    #disks
    c = type
    s = @lift( 2. *$R.*um)
    st = @lift( 2. *$R.*um/5)

    heatmap!(ax,field_x_centers.*um, field_y_centers.*um, field_C, alpha=0.2,colormap=:viridis,colorrange=(0,1))
    Colorbar(f[1,2], limits = (0, 1), label="Concentration c",colormap=:viridis)

    #vlines!(ax,field_x_centers, color="black", alpha=0.2)
    #hlines!(ax,field_y_centers, color="black", alpha=0.2)
    scatter!(ax,xt,yt, color=:blue, alpha=0.1)
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1, colormap=Reverse(:seismic))

    



    #directors
    cp = @lift( angle.($px+1im*$py) )

    


    arrows2d!(ax, x,y, px,py , color=cp,  colormap=:hsv,colorrange=(-pi,pi))

    #velocity vectors
    cv = @lift( angle.($vx+1im*$vy) )

    vp_x = @lift( cos.($cp-$cv))

    vp_y = @lift( sin.($cp-$cv) )

    #arrows!(ax, x,y, vx,vy , color=cv,  colormap=:hsv,colorrange=(-pi,pi))



    display(f)
    save_path = mkpath( joinpath(save_folder,"movie_traj.mp4"))
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

        px[] = frames[stri]["px"]
        py[] = frames[stri]["py"]

        field_C[] = frames[stri]["field_C"]

        
    end
    close(raw_data_file)
end

mainpath = joinpath(homedir(),"surfdrive","ActivePolygonClusters","simulations")#,"c_align_dp_distr","phi_0p3")

tree = construct_folder_tree_param_param_seed(mainpath)

for (param1, subdict) in tree

    for (param2, seeddict) in subdict

        if param1=="c_no_align_ndp_distr" || param1=="c_no_align_dp_distr"

        print(param1)
        print(param2)
        raw_data_file=jldopen(joinpath(mainpath,param1, param2, "simdata","raw_data.jld2"))
        make_movie(raw_data_file,joinpath(mainpath,param1, param2,"movies"))
        end

    end
end
#base_folder = "/Users/kammeraat/surfdrive/ActivePolygonClusters/simulations/align_phoretic/phi_0p1/"

base_folder = "/Users/kammeraat/surfdrive/ActivePolygonClusters/simulations/for_inference_v9/phi_0p01/"

raw_data_file = jldopen(joinpath(base_folder,"simdata","raw_data.h5"),"r",iotype=IOStream )
make_movie(raw_data_file,joinpath(base_folder,"movies"))