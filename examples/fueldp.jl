include("../src/Engine.jl")
using Random, Distributions

function simulation()

    external_forces = [ABP_2d_angular_noise()]

    pair_forces = [soft_disk_force()]

    field_forces = [field_propulsion_force(2e-2,0.01)]



    
    
    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N=1000
    Lx=150.
    Ly=150
    L=min(Lx,Ly)
    poly = 0.0000002


    initial_state = [ PolarParticle2d(i,1,0.2,0.001,[rand(Uniform(-L/2, L/2)) ,rand(Uniform(-L/2,L/2))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],1,rand(Uniform(1-poly, 1+poly)),1.,[0.,0.],[0.,0.]) for i=1:N];


    size = [Lx,Ly];

    lbin = 1
    print(lbin)

    x_bin_centers = range(start=-Lx/2, stop=Lx/2, step=lbin).+lbin/2
    y_bin_centers = range(start=-Ly/2, stop=Ly/2, step=lbin).+lbin/2

    print(x_bin_centers)

    v00=0.3

    C = ones(length(y_bin_centers), length(x_bin_centers))*v00

    field = FuelField2d(x_bin_centers,y_bin_centers,C, C.*0, C.*0)
    field_updaters = []

    system = System(size, initial_state,field,  external_forces, pair_forces, field_forces,dofevolvers, true, 1e1);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    plot_functions = (plot_sized_points!, plot_directors!, plot_velocity_vectors!)
    plot_field_functions=[plot_field_magnitude!]
    states = Euler_integrator(system, 0.1, 100000, 100000, 10, plot_functions,plot_field_functions);
    return states

end

states = simulation()