include(joinpath("..","..","src","Engine.jl"))

function simulation()

    

    pair_forces = [soft_disk_force([1,2],[1. 1. ; 1. 1.])]

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver([1,2]),overdamped_pq_evolver([1,2]))
    global_dofevolvers = []
    field_dofevolvers = []

    N=2000
    ϕ = 0.9
    poly=15e-2
    Rs = rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  sqrt(pi *sum(Rs.^2) / ϕ)
    Nal = 900
    types = vcat(ones(Int64,Nal), ones(Int64,N-Nal) .*2)

    initial_state = PolarParticle3d[ PolarParticle3d([i],[types[i]], [1], [1], [Rs[i]], [0.1], [0.01], [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];


    display(L)
    sizes = [L,L,4];
    print(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    external_forces = (ABP_3d_propulsion_force([1,2]), self_align_with_v_unit_force(1,1.0), self_align_with_v_unit_force(2,0.05),ABP_perpendicular_angular_noise([1,2],[0,0,1]))

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,2.5*(1+poly));

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.05, 1e5, Tplot=10,fps=120,plot_functions=(plot_disks_orientation!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim;

end


sim = simulation()  

@profview sim = simulation() 

@profview_allocs sim = simulation()  

 