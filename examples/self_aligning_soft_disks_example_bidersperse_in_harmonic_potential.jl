include(joinpath("..","src","Engine.jl"))

function simulation()

    external_forces = (ABP_3d_propulsion_force(1), self_align_with_v_force(1,1),ABP_perpendicular_angular_noise(1,[0,0,1]),external_harmonic_force(1,2e-3))

    pair_forces = [soft_disk_force(1,1)]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N1=500
    N2 = 500
    ϕ = 0.1
    L =  sqrt(N2 *  π * 1^2 / ϕ)
    initial_state = PolarParticle3d[PolarParticle3d(i,1, 1, 1, 1, 0.1, 0.001, [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N1 ];
    for i in N1+1:N1+N2
        push!(initial_state,PolarParticle3d(i,1, 1, 1, 2, 0.3, 0.001, [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) )

    end
    
    size = [L,L,3.];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,5.);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,5e-2, 1e5, Tplot=2e1,fps=120, plot_functions=(plot_sized_points!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end


sim = simulation()

