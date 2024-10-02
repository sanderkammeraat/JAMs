include("Engine.jl")
include("Forces.jl")
include("Particles.jl")
include("DOFEvolvers.jl")
include("LivePlottingFunctions.jl")
using Random, Distributions



dofevolvers = [overdamped_x_evolver!,overdamped_v_evolver!, overdamped_f_evolver!, overdamped_θ_evolver!, overdamped_ω_evolver!]


#Initialize state
N=200
forces = []

push!(forces, Force("ABP propulsion", "external", Dict(),contribute_2d_ABP_propulsion_force!))
#push!(forces, Force("ABP angular noise", "external", Dict(), contribute_2d_ABP_angular_noise!))

push!(forces, Force("Translatational Swarmalator force", "pair", Dict("N"=>sqrt(N), "J"=>1.), contribute_swarm_pos_force! ))
push!(forces, Force("Angular Swarmalator force", "pair", Dict("N"=>sqrt(N), "K"=>-0.1), contribute_swarm_angular_force! ))
L=100.
initial_state = [ Swarmalator(i,0.01,0,[rand(Uniform(L/3, 2*L/3)) ,rand(Uniform(L/3, 2*L/3))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],1.0,[0.,0.],[0.,0.]) for i=1:N];


size = [L,L]

system = System(size, initial_state, forces, dofevolvers, false);


#%%
states = Euler_integrator(system, 0.1, 1000, 10000, 10, plot_Swarmalators!);

