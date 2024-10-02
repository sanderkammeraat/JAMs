include("Engine.jl")
include("Forces.jl")
include("Particles.jl")
include("DOFevolvers.jl")
include("LivePlottingFunctions.jl")

using Random, Distributions


forces = []

push!(forces, Force("ABP propulsion", "external", Dict(),contribute_2d_ABP_propulsion_force!))
push!(forces, Force("ABP angular noise", "external", Dict(), contribute_2d_ABP_angular_noise!))
push!(forces, Force("soft disk repulsion", "pair",Dict(), contribute_soft_disk_force!))

dofevolvers = [overdamped_x_evolver!,overdamped_v_evolver!, overdamped_f_evolver!, overdamped_θ_evolver!, overdamped_ω_evolver!]


#Initialize state
N=100
L=20.
initial_state = [ PolarParticle2d(i,1,0.3,0.01,[rand(Uniform(0, L)) ,rand(Uniform(0, L))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],1.0,1.,1.,[0.,0.],[0.,0.]) for i=1:N];


size = [L,L];

system = System(size, initial_state, forces, dofevolvers, true);


#Run integration
states = Euler_integrator(system, 0.1, 1000, 1000, 10, plot_disks!);


