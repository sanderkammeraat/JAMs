include(joinpath("..","src","Engine.jl"))
function simulation()

    external_forces = ( ABP_3d_propulsion_force(1), self_align_with_v_force(1,0.5),ABP_perpendicular_angular_noise(1,[0,0,1]))

    

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=100
    ϕ = 0.2
    r=1.
    R =  sqrt(N * r^2 / ϕ)
    initial_state = [ PolarParticle3d([i], [1], [1], [1], [r], [0.1], [0.08], [rand(Uniform(-2*R/3, 2*R/3)) , rand(Uniform(-2*R/3,2*R/3)),0],[0,0,0],[0.,0.,0.], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];
    Nb = floor(2*pi*R/(2*r))
    print(Nb)
    j=1
    for i in N+1:N+1+Nb-1

        push!(initial_state,PolarParticle3d([i],[2], [1],[1], [r], [0], [0.01], [R*cos(2*pi/Nb*j) , R*sin(2*pi/Nb*j),0],[0,0,0],[0.,0.,0.], [0,0,0],[0,0,0],push!(normalize(rand(Normal(0),2)),0).*[1,1,0],[0,0,0],[0,0,0]))
        j+=1
    end

    pair_forces = (soft_disk_force([1, 2],[1 1; 1 0]),periodic_chain_force(2,1,0,N+1,N+1+Nb-1))

    size = [5*R+2*r,5*R+2*r,1.];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,r*4);

    #Run integration
    sim = Euler_integrator(system,1e-2, 1e5, Tplot=5e1, fps=120, plot_functions=(plot_sized_points!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end

sim = simulation()

    





