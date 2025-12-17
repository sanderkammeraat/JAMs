
include(joinpath("..","src","Engine.jl"))
include("AnalysisFunctions.jl")
include("AnalysisPipeline.jl")
using GLMakie
GLMakie.activate!()


# rp, sp, ap = auto_analysis_dir(base_folder, "sa_raw_data.h5"; support_raw_data_file_name_pattern = "ra_raw_data.h5")

# custom_analysis_function = run_sa_analysis!

# run_multithreaded_analysis(rp, ap,custom_analysis_function, support_raw_data_file_paths=sp)

base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/free/phi_1.0/N_2000/"
rp, mp =  auto_movie_dir(base_folder, "sa_raw_data.h5")



@views function make_movie_sa(raw_data_file,movie_save_path)
    save_tax = raw_data_file["integration_info"]["save_tax"]

    J = raw_data_file["system"]["forces"]["external"]["self_align_with_v_unit_force"]["β"]

    
    frame_numbers = 1:10:length(save_tax)

    frames = raw_data_file["frames"]

    Nint = length(extract_frame_data_for_type("id",1,frames["1"]))
    Dr = frames["1"]["Dr"][1]
    v0 = frames["1"]["v0"][1]   
    t = Observable(0.)


    x = Observable(frames["1"]["xuw"])
    y = Observable(frames["1"]["yuw"]) 

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
    f = Figure(size=(1500,1500));
    ax = Axis(f[1,1], aspect=DataAspect(),title = @lift("t = $(round($t, digits = 1))"), xlabel="x", ylabel="y");

    #ax2 = Axis(f[1,2],title = @lift("t = $(round($t, digits = 1)), Dr=$Dr, J=$J "), xlabel="t", ylabel="std(px), mean(px)");


    #disks
    c = @lift( angle.($px+1im*$py) )
    s = @lift( 2. *$R)

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
    #plotpx = @lift(scaleup.*$(px)[type.==1]) 
    #plotpy = @lift(scaleup.*$(py)[type.==1]) 
    #scatter!(ax,plotpx,plotpy, color=cp,  colormap=:hsv,colorrange=(-pi,pi))
    #lines!(ax2, t_stdp)
    #lines!(ax2, t_mp)
    #colsize!(f.layout, 1, Relative(3/4))
    #ylims!(ax2, (-1,1))

    k =raw_data_file["system"]["forces"]["pair"]["soft_disk_force"]["karray"]

    #Interior particles
    Nint = sum(type .== 1)
    ϕ = 1.3
    #Check if all radii are the same, if so do, else , hardcoded 0.15, because I did not store the polydispersity of the initial conditions in the  analysis file
    poly =  all( frames["1"]["R"] .== frames["1"]["R"][1]) ? 0. : 0.15


    tag = Dict("ϕ"=>ϕ, "v0"=> v0, "Nint"=> Nint, "poly"=>poly, "k"=>k, "Dr"=>Dr, "J"=>J)
    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)
    display(f)
    record(f, movie_save_path, frame_numbers, visible=true, compression = 28) do i 

        stri = string(i)
        t[] = frames[stri]["t"]

        x[] = frames[stri]["xuw"]
        y[] = frames[stri]["yuw"]

        vx[] = frames[stri]["vx"]
        vy[] = frames[stri]["vy"]

        px[] = frames[stri]["px"]
        py[] = frames[stri]["py"]

        t_stdp[] = push!(t_stdp[], Point2f(frames[stri]["t"], std( px[][type.==1]) ) )

        t_mp[] = push!(t_mp[], Point2f(frames[stri]["t"], mean( px[][type.==1]) ) )

        #xlims!(ax2, (0,1e-3+frames[stri]["t"]))
        
    end
end


#movie_single(rp[1], mp[1], make_movie_sa)

movie_single(rp[200], "/Users/kammeraat/test_uw_movie/movie.mp4", make_movie_sa)

# run_sequential_movie(rp, mp, make_movie_sa)


