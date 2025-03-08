include(joinpath("..","src","Engine.jl"))

function simulation()

    external_forces = ( ABP_3d_propulsion_force([1,2]), self_align_with_v_force([1,2],1.),ABP_perpendicular_angular_noise([1,2],[0,0,1]))

    pair_forces = [soft_disk_force([1,2,3],[1 1 2; 1 1 2; 2 2 0])]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N1=1000
    N2=1000
    ϕ = 0.5
    r1=1.
    r2=2.
    v01= 1
    v02=2.
    R =  sqrt((N1 * r1^2+N2 * r2^2) / ϕ)
    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[ PolarParticle3d(i, 1, 1, 1, r1, v01, 0.01, [rand(Uniform(-2*R/3, 2*R/3)) , rand(Uniform(-2*R/3,2*R/3)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N1 ];
    for i in N1+1:N1+N2
        push!(initial_state,PolarParticle3d(i, 2, 1, 1, r2, v02, 0.01, [rand(Uniform(-2*R/3, 2*R/3)) , rand(Uniform(-2*R/3,2*R/3)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]))
    end
    rb = r2
    Nb = floor(2*pi*R/rb)
    j=1
    for i in N1+N2+1:N1+N2+1+Nb

        push!(initial_state,ConfinedPolarParticle3d(i,3, 1, 1, rb, 0, 0.01, [R*cos(2*pi/Nb*j) , R*sin(2*pi/Nb*j),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],push!(normalize(rand(Normal(0),2)),0).*[1,1,0],[0,0,0],[0,0,0]))
        j+=1
    end
    size = [2*R+2*r2,2*R+2*2,2*r2];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,r2*4);

    #Run integration
    sim = Euler_integrator(system,1e-1, 1e5, Tplot=5e0, fps=120,plot_functions=(plot_sized_points!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end


sim = simulation()

