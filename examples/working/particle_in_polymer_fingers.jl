include(joinpath("..","..","src","Engine.jl"))


function simulation()

    pair_forces = (soft_disk_force([1,2], [1 1 ; 1 .2]),ring_polymer_harmonic_bend_force(1,0.2), ring_polymer_harmonic_stretch_force(1,1.))

    

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver([1,2]),overdamped_pq_xyc_evolver([1,2]))
    global_dofevolvers = []
    field_dofevolvers = []
    N_in_pol=100
    
    l=2
    R=N_in_pol*l/2/pi

    pf=1

    N = pf*R^2/1^2

    L=3*R
    v0=0.3
    Dr=0.01
    J=0.0
    external_forces=(ABP_3d_propulsion_force(2), ABP_perpendicular_angular_noise(2,[0, 0 ,1]), self_align_with_v_unit_force(2,J))
    
    initial_state =  PolarPolymerParticle3d[PolarPolymerParticle3d([id],[1],[1],[id],[N_in_pol], [1], [1], [1], [0.3], [0.01], [R*cos(2*pi/N_in_pol*id) , R*sin(2*pi/N_in_pol*id),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([1,0,0]),[0,0,0],[0,0,0]) for id=1:N_in_pol];
    for i=1:N

        id = N_in_pol+i
        initial_state=push!(initial_state,PolarPolymerParticle3d([id],[2],[2],[1],[1], [1], [1], [1], [v0], [Dr], [rand(Uniform(-R/2, R/2)),rand(Uniform(-R/2, R/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]))
    end

    sizes = [1.5*L,1.5*L,2.];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,6.);
    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.01, 1e5, Tplot=1, fps=60, plot_functions=(plot_sized_points!, plot_directors!, plot_velocity_vectors!), plotdim= 2); 
    return sim

end


sim = simulation()