include("Engine.jl")
include("Forces.jl")
include("Particles.jl")
include("DOFEvolvers.jl")
using Random, Distributions

forces = []

push!(forces, Force("ABP propulsion", "external", contribute_2d_ABP_propulsion_force!))
push!(forces, Force("ABP angular noise", "external", contribute_2d_ABP_angular_noise!))

dofevolvers = [overdamped_x_evolver!,overdamped_v_evolver!, overdamped_f_evolver!, overdamped_θ_evolver!, overdamped_ω_evolver!]


#Initialize state
N=100
L=100.
initial_state = [ PolarParticle2d(i,1,0.1,0.1,[L/2 ,L/2],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],1.,[0.,0.],[0.,0.]) for i=1:N];


size = [L,L];

system = System(size, initial_state, forces, dofevolvers, true);


#%%
states = Euler_integrator(system, 0.01, 1000, 100, 100);

#%%

anim = @animate for (i, state) in pairs(states)
    x = [p_i.x[1] for p_i in state]
    y = [p_i.x[2] for p_i in state]
    c = [p_i.id for p_i in state]
    scatter(x,y, xlimits = (0,system.sizes[1]), ylimits=(0, system.sizes[2]), legend=false, zcolor=c, color=:hawaii, aspect_ratio = :equal)
end

gif(anim, "anim_fps15.gif", fps = 15)

#%%
@profview Euler_integrator(system, 0.01, 500, 100)