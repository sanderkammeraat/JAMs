include(joinpath("..","..","src","Engine.jl"))

function simulation()


    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,1.2,0.3,0.3,true), pair_nematic_alignment_force(1,2.5,0.1))
    #type, torque, rfact, kpar, kper
    #1.3 -1 0 0.3
    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,true,1.3, 1, 0., 0.3), pair_nematic_alignment_force(1,2.5,0.15))
    pair_forces = (soft_disk_force(1,1),polymer_harmonic_bend_force(1,0.1), polymer_harmonic_stretch_force(1,1),polymer_align_director_tangent_force(1,10), polymer_pairAN_force(1,false,1.3, -1, 0., 0.5))

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []

    N=2000
    ϕ = 1.0
    poly=15e-6
    Rs = 1
    display(size(Rs))

    
    Npols = 40
    N_in_pol = 20

    L =  2.5*maximum([N_in_pol,Npols])
    #PolarPolymerParticle3d id type pol_id id_in_pol pol_N
    initial_state = PolarPolymerParticle3d[];
    print(typeof(initial_state))

    id=1
    for n in 1:Npols

        for i in 1:N_in_pol

            push!(initial_state,PolarPolymerParticle3d([id],[1],[n],[i],[N_in_pol], [1], [1], [Rs], [0.3], [0.01], [-N_in_pol+2*i , -Npols  + 2*n,0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]))
            id+=1
        end
    end

    display(L)
    sizes = [L,L,4];
    print(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    #β=-1 interesting!
    external_forces =(thermal_translational_noise(1, 1*[0.001, 0.0001,0]),)

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,2*2.5*(1+poly));

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.025, 1e4, Tplot=20,fps=120,plot_functions=(plot_disks_orientation!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim;

end

sim = simulation()  

@profview sim = simulation() 

@profview_allocs sim = simulation()  

 