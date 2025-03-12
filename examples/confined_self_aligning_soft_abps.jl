include(joinpath("..","src","Engine.jl"))

function simulation()

    external_forces = ( ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,1),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = [soft_disk_force([1, 2],[1 2; 2 1])]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=1000
    ϕ = 1.2
    r=1.
    R =  sqrt(N * r^2 / ϕ)
    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[ PolarParticle3d(i, 1, 1, 1, r, 0.01, 0.0001, [rand(Uniform(-2*R/3, 2*R/3)) , rand(Uniform(-2*R/3,2*R/3)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];
    Nb = floor(2*pi*R/(2*r))
    j=1
    for i in N+1:N+1+Nb

        push!(initial_state,ConfinedPolarParticle3d(i,2, 1,1, r, 0, 0.01, [R*cos(2*pi/Nb*j) , R*sin(2*pi/Nb*j),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],push!(normalize(rand(Normal(0),2)),0).*[1,1,0],[0,0,0],[0,0,0]))
        j+=1
    end
    size = [2*R+2*r,2*R+2*r,1*r];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,2.5*r);

    #Run integration
    sim = Euler_integrator(system,1e-1, 1e6, Tplot=5e0, fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end

sim = simulation();


