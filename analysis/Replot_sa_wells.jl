
include("AnalysisPipeline.jl")

function make_movie(raw_data_file,movie_save_path)
    save_tax = raw_data_file["integration_info"]["save_tax"]

    frame_numbers = 1:1:length(save_tax)

    frames = raw_data_file["frames"]
    J = raw_data_file["system"]["forces"]["external"]["self_align_with_v_unit_force"]["β"]
    Dr = frames["1"]["Dr"][1]
    v0 = frames["1"]["v0"][1]   
    t = Observable(0.)

    field_C = Observable(frames["1"]["field_C"])
    field_x_centers = frames["1"]["field_bin_centers_x"]
    field_y_centers = frames["1"]["field_bin_centers_y"]

    points = Observable(Point2f[(frames["1"]["x"][1], frames["1"]["y"][1])])

    vx = Observable(frames["1"]["vx"])
    vy = Observable(frames["1"]["vy"]) 

    x = Observable(frames["1"]["vx"])
    y = Observable(frames["1"]["vy"]) 

    px = Observable(frames["1"]["px"])
    py = Observable(frames["1"]["py"]) 
    type = frames["1"]["type"]

    R = Observable(frames["1"]["R"]) 
    
    #Setup figure
    f = Figure(size=(1000,1000))
    ax = Axis(f[1,1], aspect=DataAspect(),title = @lift("t = $(round($t, digits = 1))"), xlabel="x", ylabel="y");
    xlims!(ax, low=-2.5, high=2.5)
    ylims!(ax, low=-1.5, high=1.5)

    #potential 
    Cmax = maximum(field_C[])
    Cmin = minimum(field_C[])
    surface!(ax,field_x_centers,field_y_centers,field_C,colormap=Reverse(:gist_rainbow), alpha=0.5, colorrange=(Cmin, Cmax), color=field_C,shading=NoShading)
    Colorbar(f[1,2], limits = (Cmin, Cmax), label="Pot. energy U",colormap=Reverse(:gist_rainbow))

    scatter!(ax,points, color="black", alpha=0.2)
    meshscatter!(ax,x,y)

    #directors
    cp = @lift( angle.($px+1im*$py) )

    #arrows!(ax, x,y, px,py , color=cp,  colormap=:hsv,colorrange=(-pi,pi))

    #velocity vectors
    cv = @lift( angle.($vx+1im*$vy) )


    arrows!(ax, x,y,0, vx,vy,0 , color=cv,  colormap=:hsv,colorrange=(-pi,pi))
    tag = Dict("v0"=> v0, "Dr"=>Dr, "J"=>J)
    Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)


    display(f)
    record(f, joinpath(movie_save_path,), frame_numbers; visible=true, framerate=60) do i 

        stri = string(i)
        t[] = frames[stri]["t"]
        #tplot[]=push!(tplot[], frames[stri]["t"] )

        new_point = Point2f(frames[stri]["x"][1], frames[stri]["y"][1])
        

        points[] = push!(points[], new_point)

        x[] = frames[stri]["x"]
        y[] = frames[stri]["y"]

        vx[] = frames[stri]["vx"]
        vy[] = frames[stri]["vy"]

        # px[] = frames[stri]["px"]
        # py[] = frames[stri]["py"]

        field_C[] = frames[stri]["field_C"]

        
    end
    display(movie_save_path)
    close(raw_data_file)
end

base_folder = "/Users/kammeraat/dwsa/single/simdata/"

rp, mp =  auto_movie_dir(base_folder, "raw_data.h5")

todo =  occursin.("v0_0.4/", rp)

rpt = rp#[todo]
mpt = mp#[todo]

run_sequential_movie(rpt, mpt, make_movie, only_seed=nothing)