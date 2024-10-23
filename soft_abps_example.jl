include("Engine.jl")
include("LivePlottingFunctions.jl")

using Random, Distributions

function simulation()

    external_forces = (ABP_2d_propulsion_force(), ABP_2d_angular_noise())

    pair_forces = [soft_disk_force()]

    
    
    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N=500
    ϕ = 0.8
    L=sqrt(N*pi/ϕ)

    poly = 0.2


    initial_state = [ PolarParticle2d(i,1,0.2,0.01,[rand(Uniform(0, L)) ,rand(Uniform(0, L))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],3,rand(Uniform(1-poly, 1+poly)),1.,[0.,0.],[0.,0.]) for i=1:N];


    size = [L,L];

    system = System(size, initial_state, external_forces, pair_forces , dofevolvers, true);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    states = Euler_integrator(system, 0.1, 1000, 100000, 10, plot_polar_disks!);
    return states

end

states = simulation()


