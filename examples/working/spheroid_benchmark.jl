include(joinpath("..","..","src","Engine.jl"))

function random_outside_disk(N,R,L)
    x = zeros(N)
    y = zeros(N)

    for i in eachindex(x)
        outside_R = false
        while ~outside_R

            xi = rand(Uniform(-L/2, L/2))
            yi = rand(Uniform(-L/2, L/2))
            outside_R = xi^2 + yi^2 > R^2
            
            if outside_R
                x[i] = xi 
                y[i] = yi
            end
        end
    end
    return x,y
end

function random_inside_disk(N,R,L)
    x = zeros(N)
    y = zeros(N)
    for i in eachindex(x)
        inside_R = false
        while ~inside_R

            xi = rand(Uniform(-L/2, L/2))
            yi = rand(Uniform(-L/2, L/2))
            inside_R = xi^2 + yi^2 < R^2
            
            if inside_R
                x[i] = xi 
                y[i] = yi
            end
        end
    end
    return x,y
end

function injection()

    N_ECM=2000
    ϕ_ECM = 0.5
    
    poly=0.15
    #Size of particles
    Rs_ECM = rand(Uniform(1-poly, 1+poly),N_ECM)
    

    N_sph=600
    ϕ_sph = 0.9
    Rs_sph = rand(Uniform(1-poly, 1+poly),N_sph)


    Rdisk = sqrt(sum(Rs_sph.^2) / ϕ_sph)

    L =  sqrt(pi *sum(Rs_ECM.^2) / ϕ_ECM + pi*Rdisk^2)

    
    x_ECM, y_ECM = random_outside_disk(N_ECM,Rdisk,L)

    x_sph, y_sph = random_inside_disk(N_sph,Rdisk,L)

    initial_state = PolarParticle3d[ PolarParticle3d([i],[2], [1], [1], [Rs_ECM[i]], [0.0], [0.00], [x_ECM[i] , y_ECM[i],0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N_ECM ]

    for i in 1: N_sph
        push!(initial_state,PolarParticle3d([N_ECM+i],[1], [1], [1], [Rs_sph[i]], [0.1], [0.01], [x_sph[i] , y_sph[i],0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]))
    end

    external_forces = []

    pair_forces = (soft_atre_type_force([1,2],[1 0.4 ; 0.4 1 ],[1 1 ; 1 1 ]*0.15),)

    local_dofevolvers = (overdamped_pq_evolver([1,2]),overdamped_xvf_evolver([1,2]))
    global_dofevolvers = []
    sizes = [120,120,5];
    initial_field_state=[]
    field_forces = []
    field_updaters = []
    field_dofevovers = []
    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers, field_dofevovers, true,2.5);

    #Run integration
    sim = Euler_integrator(system,1e-1, 5000*1e-1, Tplot=nothing, fps=120, plot_functions=(plot_disks_type!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim
end

inj = injection()


function simulation(inj)

    external_forces = (ABP_perpendicular_angular_noise([1],[0,0,1]),)
    pair_forces = (soft_atre_type_force([1,2],[1 0.4 ; 0.4 1 ],[1 1 ; 1 1 ]*0.15), pairABP_force([1,2],1.5,[1 1; 1 1],false), pair_polar_alignment_force(1, 2.5, 0.1))

    local_dofevolvers = (overdamped_pq_evolver([1,2]),)
    global_dofevolvers =(overdamped_pairdis_evolver(1,1.),)
    field_dofevolvers = []

    sizes = inj.system.sizes
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    initial_particle_state =deepcopy(inj.final_particle_state)

    #Modify initial state
    for (i, p_i) in pairs(initial_particle_state)
        #reinitialize forces 
        p_i.f.*=0
        p_i.q.*=0

        if p_i.type[1]==1
            p_i.Dr[1] = 0.01/2
            p_i.v0[1] = 0.5
        end

    end

    system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers,field_dofevolvers,true,3.);

    #Run integration # 100000*1e-2
    sim = Euler_integrator(system,1e-2,100000*1e-2, Tplot=10, fps=120,Tsave=nothing, plot_functions=(plot_disks_type!, plot_velocity_vectors!,plot_directors!), plotdim=2, save_folder_path=joinpath(homedir(),"sph","benchmark"),save_functions=[save_2d_polar_p!]); 
    return sim
end

display(Threads.nthreads())
simulation(inj)

@profview simulation(inj)

