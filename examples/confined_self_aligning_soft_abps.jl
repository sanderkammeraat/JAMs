include("../src/Engine.jl")

function simulation()

    external_forces = ( ABP_3d_propulsion_force(), self_align_with_v_force(0.5),ABP_perpendicular_angular_noise([0,0,1]))

    pair_forces = [ new_soft_disk_force(2)]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=1000
    ϕ = 1.2
    r=1.
    R =  sqrt(N * r^2 / ϕ)
    initial_state = Union{NewPolarParticle3d,ConfinedPolarParticle3d}[ NewPolarParticle3d(i, 1, 1, r, 0.3, 0.08, [rand(Uniform(-2*R/3, 2*R/3)) , rand(Uniform(-2*R/3,2*R/3)),0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0]) for i=1:N ];
    Nb = floor(2*pi*R/r)
    j=1
    for i in N+1:N+1+Nb

        push!(initial_state,ConfinedPolarParticle3d(i, 1, 1, r, 0, 0.01, [R*cos(2*pi/Nb*j) , R*sin(2*pi/Nb*j),0],[0,0,0], [0,0,0],[0,0,0],push!(normalize(rand(Normal(0),2)),0).*[1,1,0],[0,0,0]))
        j+=1
    end
    size = [2*R+2*r,2*R+r,2*R+2*r];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,r*10);

    #Run integration
    states = Euler_integrator(system,1e-1, 1e5, 1e10, 5e0, 120,(new_plot_sized_points!,plot_directors!, plot_velocity_vectors!), true); 
    return states

end


states = simulation()

