include("../src/Engine.jl")


function simulation()

    external_forces = [external_harmonic_force(1e-1), ABP_3d_propulsion_force(), self_align_with_v_force(1), external_friction_force()]

    pair_forces = []

    dofevolvers = [inertial_evolver!]
s
    initial_state = [ Hexbug(1, 1, 1, 0, 0.1, 0.001, [1,1,0],[0,0,0], [0,0,0],[0,0,0],[1,0,0],[0,0,0]) ];

    L = 5.
    size = [L,L,L];

    system = System(size, initial_state, external_forces, pair_forces , dofevolvers, true);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    states = Euler_integrator(system, 1e-3, 100000, 100000, 1e3, (plot_points!, plot_directors!, plot_velocity_vectors!), true);
    return states

end


states = simulation()
