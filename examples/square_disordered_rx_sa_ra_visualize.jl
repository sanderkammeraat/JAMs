
include(joinpath("..","src","Engine.jl"))

function relaxation_step(save_folder_path; Tsave=100, Tplot=nothing)

    external_forces = []#[thermal_translational_noise(1, 0 .*[1.,1.,0])]

    pair_forces = [soft_disk_force([1, 2],[1. 1.; 1. 1.])]
    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []


   


    #First  fixed outline
    Nlin=20

    xs = range(-(Nlin-1),(Nlin-1), step=2)

    ys = range(-(Nlin-1),(Nlin-1), step=2)



    xsq = []
    ysq = []
    #Loop clockwise
    for xi in xs
        push!(xsq, xi)
        push!(ysq, ys[end])
    end

    for yi in reverse(ys)[2:end]
        push!(xsq, xs[end])
        push!(ysq, yi)
    end

    for xi in reverse(xs)[2:end]
        push!(xsq, xi)
        push!(ysq, ys[1])
    end

    for yi in ys[2:end-1]
        push!(xsq, xs[1])
        push!(ysq, yi)
    end
    print(xs)
    print(xsq)



    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[]

    r=1.0
    ϕ=0.9
    typess = []

    x = xsq
    y = ysq
    Nfix = length(x)

    print(Nfix)
    poly=0.15
    Rsfix = rand(Uniform((1-poly)*r, (1+poly)*r),Nfix)
    while mean(Rsfix)<1 || mean(Rsfix)>1+ 1e-2
        Rsfix = rand(Uniform((1-poly)*r, (1+poly)*r),Nfix)
        println(mean(Rsfix))
    end

    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[]
    id=1

    for i=1:Nfix
        push!(initial_state,ConfinedPolarParticle3d([id],[2], [1],[1], [Rsfix[i]], [0.], [0.], [x[i] , y[i],0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]))
        id+=1
    end

    N = round(Int64,( (2*Nlin)^2 - pi * sum(Rsfix.^2)/2 )*ϕ/pi )
    print(N)

    Rs = rand(Uniform((1-poly)*r, (1+poly)*r),N)
    while mean(Rs)<1 || mean(Rs)>1+ 1e-2
        Rs = rand(Uniform((1-poly)*r, (1+poly)*r),N)
        println(mean(Rs))
    end
    #Then loop over normal particles
    for i=1:N

            
        push!(initial_state, PolarParticle3d([id], [1], [1], [1], [Rs[i]], [0.], [0.], vcat(rand(Uniform(-Nlin//2,Nlin//2),2),[0.]),[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]))
        id+=1
    end


    size = [Nlin*2+3, Nlin*2+3,1];
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers, field_dofevolvers, false,2.5*r*(1+poly));

    #Run integration
    sim = Euler_integrator(system,1e-1, 2e3,Tsave=Tsave, Tplot=Tplot, fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end


function self_aligning_step(rx_step,J,v0, Dr, seed,save_folder_path; Tsave=100, Tplot=nothing)

    external_forces =( ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,J),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = (soft_disk_force([1, 2],[1 1.; 1. 1]),)

    local_dofevolvers = [overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1)]
    global_dofevolvers = []
    field_dofevolvers = []

    sizes = 2 .* rx_step.system.sizes
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    initial_particle_state =deepcopy(rx_step.final_particle_state)

    #Modify initial state
    for (i, p_i) in pairs(initial_particle_state)
        #reinitialize forces 
        p_i.f.*=0
        p_i.q.*=0
        p_i.Dr[1] = Dr
        p_i.v0[1] = v0

    end
    system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers,field_dofevolvers,false,rx_step.system.rcut_pair_global);

    #Run integration
    sim = Euler_integrator(system,1e-2, 1e5,Tsave=Tsave,seed=seed, Tplot=Tplot, fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim
end


function relax_again_step(sa_step, save_folder_path; Tsave=100, Tplot=nothing)

    external_forces =[] # [thermal_translational_noise(1, 0 .*[1.,1.,0])]

    pair_forces = [soft_disk_force([1, 2],[1. 1.; 1. 1.])]

    local_dofevolvers = [overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1)]
    global_dofevolvers = []
    field_dofevolvers = []

    sizes = sa_step.system.sizes
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    initial_particle_state =deepcopy(sa_step.final_particle_state)

    #Modify initial state
    for (i, p_i) in pairs(initial_particle_state)
        #reinitialize forces 
        p_i.f.*=0
        p_i.q.*=0
        p_i.Dr[1] = 0.
        p_i.v0[1] = 0.

    end
    system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers,field_dofevolvers,false,sa_step.system.rcut_pair_global);

    #Run integration
    sim = Euler_integrator(system,1e-1, 5e3,Tsave=Tsave,seed=nothing, Tplot=Tplot)# , fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim
end

rx_result= relaxation_step("",Tsave=nothing, Tplot=100)
sa_result=self_aligning_step(rx_result,1,0.001, 0.00,1, ""; Tsave=nothing, Tplot=100);

ra_result=relax_again_step(sa_result, ""; Tsave=nothing, Tplot=100);



