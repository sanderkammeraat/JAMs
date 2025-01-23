include("../src/Engine.jl")

function simulation()

    external_forces = [ ABP_3d_propulsion_force(), self_align_with_v_force(1),ABP_perpendicular_angular_noise([0,0,1]),external_harmonic_force(1e-3)]

    pair_forces = [ new_soft_disk_force(1)]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N1=500
    N2 = 500
    ϕ = 0.1
    L =  sqrt(N2 *  π * 1^2 / ϕ)
    initial_state = [ NewPolarParticle3d(i, 1, 1, 1, 0.2, 0.001, [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0]) for i=1:N1 ];
    for i in N1+1:N1+N2
        push!(initial_state,NewPolarParticle3d(i, 1, 1, 1, 0.2, 0.001, [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0]) )

    end
    
    size = [L,L,L];

    system = System(size, initial_state, external_forces, pair_forces , dofevolvers, true);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    states = Euler_integrator(system,1e-1, 1e5, 1e10, 1e1, (new_plot_sized_points!, plot_directors!, plot_velocity_vectors!), true); 
    return states

end


states = simulation()

