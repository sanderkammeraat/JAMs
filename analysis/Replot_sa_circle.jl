
include(joinpath("..","src","Engine.jl"))
using GLMakie




base_folder = joinpath(homedir(), "sa", "production","phi_1","tstop_2e3","simdata")
begin
begin
#load relaxation
relax_file = jldopen( joinpath(base_folder,"relaxation_step" ,"raw_data.jld2"), "r")

#Setup observables
save_tax = relax_file["integration_info"]["save_tax"]

frame_numbers = 1:length(save_tax)

frames = relax_file["frames"]


t = Observable(0.)

x = Observable(frames["1"]["x"])
y = Observable(frames["1"]["y"]) 

vx = Observable(frames["1"]["vx"])
vy = Observable(frames["1"]["vy"]) 

px = Observable(frames["1"]["px"])
py = Observable(frames["1"]["py"]) 

R = Observable(frames["1"]["R"]) 

#Setup figure
f = Figure(size=(1000,1000));
ax = Axis(f[1,1], aspect=1,title = @lift("t = $(round($t, digits = 1))"));

c = @lift( angle.($px+1im*$py) )
s = @lift( 2. *$R )
scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))
display(f)
end

record(f, joinpath(base_folder,"relaxation_step.mp4"), frame_numbers; framerate=60, visible=true) do i 

    stri = string(i)
    t[] = frames[stri]["t"]

    x[] = frames[stri]["x"]
    y[] = frames[stri]["y"]

    vx[] = frames[stri]["vx"]
    vy[] = frames[stri]["vy"]

    px[] = frames[stri]["px"]
    py[] = frames[stri]["py"]
end


#load relaxation



begin
sa_file =  jldopen( joinpath(base_folder,"self_alignment_step" ,"raw_data.jld2"), "r") 
#Setup observables
save_tax = sa_file["integration_info"]["save_tax"]

frame_numbers = 1:length(save_tax)

frames = sa_file["frames"]


t = Observable(0.)

x = Observable(frames["1"]["x"])
y = Observable(frames["1"]["y"]) 

vx = Observable(frames["1"]["vx"])
vy = Observable(frames["1"]["vy"]) 

px = Observable(frames["1"]["px"])
py = Observable(frames["1"]["py"]) 

R = Observable(frames["1"]["R"]) 

#Setup figure
f = Figure(size=(1000,1000));
ax = Axis(f[1,1], aspect=1,title = @lift("t = $(round($t, digits = 1))"));

c = @lift( angle.($px+1im*$py) )
s = @lift( 2. *$R )
scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))
display(f)
end

record(f, joinpath(base_folder,"self_alignment_step.mp4"), frame_numbers; framerate=60, visible=true) do i 

    stri = string(i)
    t[] = frames[stri]["t"]

    x[] = frames[stri]["x"]
    y[] = frames[stri]["y"]

    vx[] = frames[stri]["vx"]
    vy[] = frames[stri]["vy"]

    px[] = frames[stri]["px"]
    py[] = frames[stri]["py"]
end


begin
#load relax again step

relax_again_file = jldopen( joinpath(base_folder,"relax_again_step" ,"raw_data.jld2"), "r") 

#Setup observables
save_tax = relax_again_file["integration_info"]["save_tax"]

frame_numbers = 1:length(save_tax)

frames = relax_again_file["frames"]


t = Observable(0.)

x = Observable(frames["1"]["x"])
y = Observable(frames["1"]["y"]) 

vx = Observable(frames["1"]["vx"])
vy = Observable(frames["1"]["vy"]) 

px = Observable(frames["1"]["px"])
py = Observable(frames["1"]["py"]) 

R = Observable(frames["1"]["R"]) 

#Setup figure
f = Figure(size=(1000,1000));
ax = Axis(f[1,1], aspect=1,title = @lift("t = $(round($t, digits = 1))"));

c = @lift( angle.($px+1im*$py) )
s = @lift( 2. *$R )
scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))
display(f)
end

record(f,joinpath(base_folder,"relax_again_step.mp4"), frame_numbers; framerate=60, visible=true) do i 

    stri = string(i)
    t[] = frames[stri]["t"]

    x[] = frames[stri]["x"]
    y[] = frames[stri]["y"]

    vx[] = frames[stri]["vx"]
    vy[] = frames[stri]["vy"]

    px[] = frames[stri]["px"]
    py[] = frames[stri]["py"]
end
end

