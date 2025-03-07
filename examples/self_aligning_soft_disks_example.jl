include("../src/Engine.jl")


function simulation()

    external_forces = (ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,0.1),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = [soft_disk_force(1,1)]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=5000
    ϕ = 1.
    poly=1e-4
    Rs = rand(Uniform(1-poly, 1+poly),N)

    display(size(Rs))

    L =  sqrt(pi *sum(Rs.^2) / ϕ)
    initial_state = PolarParticle3d[ PolarParticle3d(i,1, 1, 1, Rs[i], 0.1, 0.1, [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];

    
    sizes = [L,L,4];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,2.5);

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,1e-1, 1e5, Tplot=5e0, fps=120, plot_functions=(plot_disks!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end

sim = simulation()