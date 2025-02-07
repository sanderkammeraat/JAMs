include("../src/Engine.jl")

using Random, Distributions
using JET
#using BenchmarkTools
function simulation()

    external_forces =  (ABP_2d_propulsion_force(), ABP_2d_angular_noise())

    pair_forces = [soft_disk_force()]

    
    
    dofevolvers = [overdamped_evolver!]

    #Initialize state
    N=1000
    L=100.
    initial_particle_state = [ PolarParticle2d(i,1,0.3,0.01,[rand(Uniform(0, L)) ,rand(Uniform(0, L))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],1.0,1.,1.,[0.,0.],[0.,0.]) for i=1:N];
    initial_field_state = []
    field_forces = []
    field_updaters = []

    size =  [L,L];

    system =System(size, initial_particle_state,initial_field_state,external_forces, pair_forces, field_forces, field_updaters,dofevolvers, true, 1e1);

    #Run integration
    states = Euler_integrator(system, 0.1, 1, 1000000, 10,120, [plot_points!]);
    0

end

simulation()

@time simulation()

@profview simulation()

@report_opt simulation()

