include("../src/Engine.jl")
using Random, Distributions

function simulation()

    external_forces = (ABP_3d_propulsion_force(), ABP_3d_angular_noise())

    pair_forces = [soft_disk_force()]

    
    
    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N=500
    L=20.

    poly = 0.2

    
    initial_state = [ PolarParticle3dN(i,rand(Uniform(1-poly, 1+poly)),1,0.1,0.01, rand(Uniform(-L/2, L/2),3),[0.,0.,0.],[0.,0.,0.],normalize(rand(Normal(0, 1),3)),[0,0,0],[0.],[0,0,0],1,[0.,0.,0.],[0.,0,0.]) for i=1:N];
    
    size = [L,L,L];

    system = System(size, initial_state, external_forces, pair_forces , dofevolvers, true);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    states = Euler_integrator(system, 0.1, 1000, 100000, 10, (plot_sized_points!, plot_directors!) );
    return states

end


states = simulation();


 