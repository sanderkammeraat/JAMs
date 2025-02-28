include("../src/Engine.jl")


function initial_condition()

    external_forces = ( ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,0),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = [soft_disk_force([1, 2],[1 2; 2 1])]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=2000
    poly=1e-1
    Rs = rand(Uniform(1-poly, 1+poly),N)
    ϕ = 0.9
    R =  sqrt(sum(Rs.^2 )/ ϕ)


    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[ PolarParticle3d(i, 1, 1, 1, Rs[i], 0.1, 0.001, [rand(Uniform(-2*R/3, 2*R/3)) , rand(Uniform(-2*R/3,2*R/3)),0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];
    r=1
    Nb = floor(2*pi*R/(2*r))
    j=1
    for i in N+1:N+1+Nb

        push!(initial_state,ConfinedPolarParticle3d(i,2, 1,1, r, 0, 0.01, [R*cos(2*pi/Nb*j) , R*sin(2*pi/Nb*j),0],[0,0,0], [0,0,0],[0,0,0],push!(normalize(rand(Normal(0),2)),0).*[1,1,0],[0,0,0],[0,0,0]))
        j+=1
    end
    size = [2*R+2*r,2*R+2*r,1*r];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,2.5*r);

    #Run integration
    sim = Euler_integrator(system,1e-1, 200, 1e1*200,5e0, 120,(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), 2); 

    return sim
end

ic = initial_condition()




function simulation(ic)

    end_state = ic.particle_states[end]

    pair_forces = [soft_disk_force(1,1.)]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]

    initial_state = PolarParticle3d[]

    for p_i in end_state

        if typeof(p_i)==PolarParticle3d

            p_im = deepcopy(p_i)
            push!(initial_state, p_im)

        end


    end
    
    pins  = zeros(length(initial_state),3)
    for i=eachindex(initial_state)

        pins[i,:].= initial_state[i].x

    end

    
    sizes = ic.system.sizes
    print(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    external_forces = (ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,1),ABP_perpendicular_angular_noise(1,[0,0,1]),external_harmonic_pinning_force(1,0.1,0,pins))

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,2.5);

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,1e-1, 5e2, 5e0, 5e0, 120,(plot_disks_orientation!, plot_directors!, plot_velocity_vectors!), 2); 
    return sim;

end


sim = simulation(ic)  

make_movie(sim, "/Users/kammeraat/test_JAMS/movies/","poly_disordered_save.mp4",(plot_disks_orientation!, plot_directors!, plot_velocity_vectors!),60,2)