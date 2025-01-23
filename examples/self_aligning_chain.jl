include("../src/Engine.jl")

function simulation()

    external_forces = [ ABP_3d_propulsion_force(), self_align_with_v_unit_force(1),ABP_perpendicular_angular_noise([0,0,1])]

    pair_forces = (chain_force(1e0,2),new_soft_disk_force(1e-5))

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=20
    L=4 * N*1.
    initial_state = [ NewPolarParticle3d(1, 1, 1, 1, 0*0.3, 0.1, [0 , -1,0],[0,0,0], [0,0,0],[0,0,0],normalize([0,1,0]),[0,0,0])];

    for i in 2:N-1
        v0 =  0.2/i
        push!(initial_state, NewPolarParticle3d(i, 1, 1, 1,v0, 0.2, [0 , -i*2,0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0]))
        print(v0)
    end

    #push!(initial_state, NewPolarParticle3d(N, 1, 1, 1, 0, 0.01, [0 , -N*2,0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0]))
    
    size = [L,L,L];

    system = System(size, initial_state, external_forces, pair_forces , dofevolvers, true);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    states = Euler_integrator(system,1e-3, 1e5, 1e10, 1e3, (new_plot_sized_points!, plot_directors!, plot_velocity_vectors!), true); 
    return states

end


states = simulation()

