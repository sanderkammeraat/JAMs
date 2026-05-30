include(joinpath("..","src","Engine.jl"))


function relaxation(base_folder)

    external_forces =[]# ( ABP_3d_propulsion_force(1),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = [soft_disk_force([1, 2],[1 2; 1 2])]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]

    N = 1000
    r_avg = 1
    ϕ = 1
    R = sqrt( N*pi/( ϕ * pi * r_avg^2) )
    bdens = 1.2
    Nb = floor(Int64,bdens*2*pi*R/(2*r_avg))
    
    poly = 0.3

    r = rand(Uniform(r_avg-poly/2, r_avg+poly/2), N+Nb)

    r_max = r_avg + poly

    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[ PolarParticle3d([i], [1], [1.], [1.], [r[i]], [0.2], [0.01], [rand(Uniform(-2*R/3, 2*R/3)) , rand(Uniform(-2*R/3,2*R/3)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];
  
    j=1
    for i in N+1:N+Nb

        push!(initial_state,ConfinedPolarParticle3d([i], [2], [1.],[1.], [r[i]], [0.], [0.0], [R*cos(2*pi/Nb*j) , R*sin(2*pi/Nb*j),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],push!(normalize(rand(Normal(0),2)),0),[0,0,0],[0,0,0]))
        j+=1
    end
    size = [2*R+2*r_max, 2*R+2*r_max, 1*r_max];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,2.5*r_max);

    #Run integration
    relax = Euler_integrator(system,1e-1, 1e4, Tplot=100, fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2, Tsave = 1000, save_functions = [save_2d_polar_p!], save_folder_path = joinpath(base_folder,"relaxation_step")); 
    return relax
end

function self_alignment_step(relax,base_folder)

    external_forces = ( ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,1),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = [soft_disk_force([1, 2],[1 2; 2 1])]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]


    sizes = relax.system.sizes
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    initial_particle_state =deepcopy(relax.final_particle_state)
    print(typeof(initial_particle_state))

    #Modify initial state
    for (i, p_i) in pairs(initial_particle_state)

        #reinitialize forces 
        p_i.f.*=0
        p_i.q.*=0
        p_i.v0[1]= 0.01
        p_i.Dr[1]= 0.01

        

    end


    system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,relax.system.rcut_pair_global);

    #Run integration
    sim = Euler_integrator(system,1e-2, 2e3, Tplot=1000, fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2, Tsave = 100, save_functions = [save_2d_polar_p!], save_folder_path = joinpath(base_folder,"self_alignment_step")); 
    return sim
end

function relax_again_step(self_align,base_folder)

    external_forces = []#( ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,1),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = [soft_disk_force([1, 2],[1 2; 2 1])]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]


    sizes = self_align.system.sizes
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    initial_particle_state = deepcopy(self_align.final_particle_state)

    #Modify initial state
    for (i, p_i) in pairs(initial_particle_state)

        #reinitialize forces 
        p_i.f.*=0
        p_i.q.*=0
        p_i.v0[1]= 0.01
        p_i.Dr[1]= 0.01

    end


    system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,self_align.system.rcut_pair_global);

    #Run integration
    sim = Euler_integrator(system,1e-1, 1e4, Tplot=100, fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2, Tsave = 1000, save_functions = [save_2d_polar_p!], save_folder_path = joinpath(base_folder,"relax_again_step")); 
    return sim
end

#Run the relaxation step


base_folder = joinpath(homedir(),"sa","production","phi_1","tstop_2e3","simdata")
relax = relaxation(base_folder)


#Run the self-alignment step
self_align = self_alignment_step(relax,base_folder);

relax_again = relax_again_step(self_align,base_folder)


