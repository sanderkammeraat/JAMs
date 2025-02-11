include("../src/Engine.jl")

using Random, Distributions
using JET
#using BenchmarkTools
begin
external_forces =  (ABP_2d_propulsion_force(1), ABP_2d_angular_noise(1))

pair_forces =[soft_disk_force(1,1)]



dofevolvers = [overdamped_evolver!]

#Initialize state
N=1000
L=60.
initial_particle_state = [ PolarParticle2d(i,1,1,0.3,0.01,[rand(Uniform(0, L)) ,rand(Uniform(0, L))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],1.,1.,[0.,0.],[0.,0.],[0,0]) for i=1:N];
initial_field_state =[]
field_forces =[]
field_updaters = []

size =  [L,L];

system =System(size, initial_particle_state,initial_field_state,external_forces, pair_forces, field_forces, field_updaters,dofevolvers, true, 1e1);
end
#Run integration
function run(system)
    Euler_integrator(system, 0.1, 1, 1e4,5,120);
end

run(system);

@time run(system);

@report_opt run(system)

@profview run(system)
