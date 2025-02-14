include("../src/Engine.jl")

function simulation()

    external_forces = ( ABP_3d_propulsion_force(1), self_align_with_v_force(1,0.5),ABP_perpendicular_angular_noise(1,[0,0,1]))

    

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=1000
    ϕ = 1.2
    r=1.
    R =  sqrt(N * r^2 / ϕ)
    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[ PolarParticle3d(i, 1, 1, 1, r, 0.3, 0.08, [rand(Uniform(-2*R/3, 2*R/3)) , rand(Uniform(-2*R/3,2*R/3)),0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];
    Nb = floor(2*pi*R/r)
    j=1
    for i in N+1:N+1+Nb

        push!(initial_state,PolarParticle3d(i,2, 1,1, r, 0, 0.01, [R*cos(2*pi/Nb*j) , R*sin(2*pi/Nb*j),0],[0,0,0], [0,0,0],[0,0,0],push!(normalize(rand(Normal(0),2)),0).*[1,1,0],[0,0,0],[0,0,0]))
        j+=1
    end

    pair_forces = (soft_disk_force([1, 2],[1 2; 2 1]),periodic_chain_force(2,1,r*0,N+1,N+1+Nb))

    size = [4*R+2*r,4*R+2*r,2*r];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,r*2.5);

    #Run integration
    sim = Euler_integrator(system,1e-1, 1e5, 1e10,5e0, 120,(plot_sized_points!,plot_directors!, plot_velocity_vectors!), 2); 
    return sim

end


sim = simulation()

    





