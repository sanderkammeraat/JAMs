include("../src/Engine.jl")
using Random, Distributions

function simulation()

    external_forces = (ABP_2d_propulsion_force(1), ABP_2d_angular_noise(1))

    pair_forces = [soft_disk_force(1,1)]

    
    
    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N=1000
    ϕ = 0.3
    L=sqrt(N*pi/ϕ)
    poly = 1e-3

    initial_state = [ PolarParticle2d(i,1,1,0.3,0.0001,[rand(Uniform(-L/2, L/2)) ,rand(Uniform(-L/2, L/2))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],rand(Uniform(1-poly, 1+poly)),1.,[0.,0.],[0.,0.],[0,0]) for i=1:N];


    size = [L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,4e0);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system, 0.1,1e5, 1e8, 5,120, (plot_sized_points!, plot_directors!, plot_velocity_vectors!));
    return sim
end


sim = simulation();