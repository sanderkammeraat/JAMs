using Observables
using GLMakie

include("loaddata.jl")

GLMakie.activate!()

@views function make_movie(raw_data_file, save_folder)
    save_tax = raw_data_file["integration_info"]["save_tax"]

    sizes =  raw_data_file["system/sizes"]
    #p_value = raw_data_file["system/forces/pair/polymer_pairAN_force/parray"]
    #kpar_value = raw_data_file["system"]["forces"]["pair"]["polymer_pairAN_force"]["k_par"]
    #kperp_value = raw_data_file["system"]["forces"]["pair"]["polymer_pairAN_force"]["k_per"]
    frame_numbers = 1:1:length(save_tax)

    frames = raw_data_file["frames"]
    #Dr = frames["1"]["Dr"][1]

    t = Observable(0.)

    x_wrapped = zeros(length(frames), length(frames["1/xuw"]))
    y_wrapped = zeros(length(frames), length(frames["1/yuw"]))

    for i in 1:length(frames)
        x_wrapped[i,:] += (frames["$i/xuw"] .+ sizes[1]/2) .% sizes[1] .- sizes[1]/2
        y_wrapped[i,:] += (frames["$i/yuw"] .+ sizes[2]/2) .% sizes[2] .- sizes[2]/2
    end

    x = Observable(x_wrapped[1,:])
    y = Observable(y_wrapped[1,:]) 

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

    pxminus = @lift(.-1 * $px)
    pyminus = @lift(.-1 * $py)
    
    arrows!(ax, x, y, px, py, color=cp, colormap=:hsv, colorrange=(-pi,pi))
    arrows!(ax, x, y, pxminus, pyminus, color=cp, colormap=:hsv,colorrange=(-pi,pi))

    #velocity vectors
    cv = @lift( angle.($vx+1im*$vy) )

    vp_x = @lift( cos.($cp-$cv))

    vp_y = @lift( sin.($cp-$cv) )

    #arrows!(ax, x,y, vx,vy , color=cv,  colormap=:hsv,colorrange=(-pi,pi))



    display(f)
    save_path = mkpath(joinpath(save_folder,"sim_movie.mp4"))
    record(f,save_path, frame_numbers; visible=true) do i 

        stri = string(i)
        t[] = frames[stri]["t"]
        #tplot[]=push!(tplot[], frames[stri]["t"] )

        x[] = x_wrapped[i,:]
        y[] = y_wrapped[i,:]

        vx[] = frames[stri]["vx"]
        vy[] = frames[stri]["vy"]

        px[] = frames[stri]["px"]
        py[] = frames[stri]["py"]

        
    end
    close(raw_data_file)
end

sim_folder_name = "sim_data"


p = [0.04, 0.06, 0.08, 0.1, 0.13, 0.15, 0.2, 0.4]
kactive = [(-1., 0.), (1., 0.), (0., 1.), (0., -1.), (1/sqrt(2), 1/sqrt(2)),(-1/sqrt(2), 1/sqrt(2)),(1/sqrt(2), -1/sqrt(2)),(-1/sqrt(2), -1/sqrt(2))]



for p_value in p
for (kpar_value, kperp_value) in kactive

    parameters = "p_$p_value,kpar_$kpar_value,kperp_$kperp_value"

    #for windows
    base_folder = joinpath("E:", "martin", sim_folder_name, parameters)

    #for linux
    #base_folder = joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name, parameters)

    if !isfile(joinpath(base_folder, "movie", "sim_movie.mp4"))

        raw_data_file = load_file(joinpath(base_folder, "raw_data.h5"))
        make_movie(raw_data_file,joinpath(base_folder,"movie"))
        close(raw_data_file)
        GLMakie.closeall()
        
    end
    
end
end