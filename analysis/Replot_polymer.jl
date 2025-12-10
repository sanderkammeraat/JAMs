
include(joinpath("..","src","Engine.jl"))
include("AnalysisFunctions.jl")
include("AnalysisPipeline.jl")
using GLMakie

GLMakie.activate!()

@views function make_movie(raw_data_file,save_folder)
    save_tax = raw_data_file["integration_info"]["save_tax"]

    k_par = raw_data_file["system"]["forces"]["pair"]["polymer_pairAN_force"]["k_par"]
    k_per = raw_data_file["system"]["forces"]["pair"]["polymer_pairAN_force"]["k_per"]
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
    pol_id = Observable(frames["1"]["pol_id"]) 

    R = Observable(frames["1"]["R"]) 
    
    #Setup figure
    f = Figure(size=(1000,1000))#Figure(size=(2000,2000));
    ax = Axis(f[1,1], aspect=DataAspect(),title = @lift("t = $(round($t, digits = 1))"), xlabel="x", ylabel="y");


    #disks
    c = pol_id
    s = @lift( 2. *$R)

    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1, colormap=:rainbow1)



    #directors
    cp = @lift( angle.($px+1im*$py) )

    


    #arrows!(ax, x,y, px,py , color=cp,  colormap=:hsv,colorrange=(-pi,pi))

    #velocity vectors
    cv = @lift( angle.($vx+1im*$vy) )

    vp_x = @lift( cos.($cp-$cv))

    vp_y = @lift( sin.($cp-$cv) )

    #arrows!(ax, x,y, vx,vy , color=cv,  colormap=:hsv,colorrange=(-pi,pi))



    display(f)
    save_path = mkpath( joinpath(save_folder,"polymer_k_par_$(k_par)_k_per_$(k_per).mp4"))
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

        
    end
    close(raw_data_file)
end

base_folder= joinpath("/data1/martin/sim_data/p_0.1/")


raw_data_file=h5open(joinpath(base_folder,"raw_data.h5"),"r")

make_movie(raw_data_file,joinpath(base_folder,"movie"))
