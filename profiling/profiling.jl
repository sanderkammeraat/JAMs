include(joinpath("..","src","Engine.jl"))

begin
external_forces =  (ABP_2d_propulsion_force(1), ABP_2d_angular_noise(1))

pair_forces = [soft_disk_force(1,1)]

dofevolvers = [overdamped_evolver!]

#Initialize state
N=10000
ϕ = 0.9
L=sqrt(N*pi/ϕ)
seed = 33
Random.seed!(seed)
initial_particle_state = [ PolarParticle2d(i,1,1,0.3,0.01,[rand(Uniform(-L/2, L/2)) ,rand(Uniform(-L/2, L/2))],[0.,0.],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],1.,1.,[0.,0.],[0.,0.],[0,0]) for i=1:N];
initial_field_state =[]
field_forces =[]
field_updaters = []

sizes =  [L,L];

system =System(sizes, initial_particle_state,initial_field_state,external_forces, pair_forces, field_forces, field_updaters,dofevolvers, true, 2.5);
end

#check integration
Euler_integrator(system, 0.1,1,seed=33,Tplot=10,fps=120, plot_functions=[plot_disks!]);

#Run integration
function runsim(system)
    Euler_integrator(system, 0.01, 0.1,seed=33)
end

runsim(system);

@time runsim(system);

@profview runsim(system);

@profview_allocs runsim(system);

