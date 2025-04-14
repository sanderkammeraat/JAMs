

#Run this script with one thread from terminal, n will run the different parameter space points over n cores
using ArgParse
using Distributed
simargs = ArgParseSettings()

@add_arg_table simargs begin
    "--ncores", "-n"
        arg_type=Int64

end

parsed_params = parse_args(simargs)
n=parsed_params["ncores"]

addprocs(n)
@everywhere include(joinpath(pwd(),"src","Engine.jl"))

@everywhere function relaxation_step(save_folder_path; Tsave=100, Tplot=nothing)

    external_forces =[thermal_translational_noise(1, 0 .*[1.,1.,0])]

    pair_forces = [soft_disk_force([1, 2],[1. 1.; 1. 1.])]
    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []
    #First make stair
    Nlin=4
    Nrows = 2*Nlin
    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[]
    xs = []
    ys = []
    r=1.0
    ϕ=0.5
    l = 2* sqrt(pi*sqrt(3)/(6*ϕ))# 2*r


    typess = []

    for row in 1:Nlin

        push!(xs, [xi-(Nlin+1)*l/2-(row-1)*l/2 for xi in l.*range(1,Nlin+row-1) ])

        push!(ys, [-(Nlin-row)/2*l*sqrt(3) for n in xs[row] ])

        if row==1
            rowtypes=2*ones(Nlin)

        else
            rowtypes = ones(length(xs[row]))
            rowtypes[1]=2
            rowtypes[end]=2

        end
        push!(typess,rowtypes )

        
    end

    for row in 1:Nlin-1
        push!(xs, xs[Nlin-row])
        push!(ys, ys[Nlin][1].-ys[Nlin-row])

        push!(typess, typess[Nlin-row])

    end

    x = vcat(xs...)
    y = vcat(ys...)
    types= vcat(typess...)
    N = length(x)
    println(N)

    poly=0.15
    Rs = rand(Uniform((1-poly)*r, (1+poly)*r),N)
    while mean(Rs)<1 || mean(Rs)>1+ 1e-2
        Rs = rand(Uniform((1-poly)*r, (1+poly)*r),N)
        println(mean(Rs))
    end

    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[]
    id=1
    #First loop over normal particles
    for i=1:N

        if types[i]==1
            
            push!(initial_state, PolarParticle3d([id], [types[i]], [1], [1], [Rs[i]], [0.], [0.], vcat(rand(Uniform(-l*Nlin/3,l*Nlin/3),2),[0.]),[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]))
            id+=1
        end
    end

    for i=1:N

        if types[i]==2
            push!(initial_state,ConfinedPolarParticle3d([id],[2], [1],[1], [Rs[i]], [0.], [0.], [x[i] , y[i],0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]))
            id+=1
        end
    end

    size = [Nrows*l+2*l,Nrows*l+2*l,1];
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers, field_dofevolvers, false,2.5*r*(1+poly));

    #Run integration
    sim = Euler_integrator(system,1e-1, 5e4,Tsave=Tsave, Tplot=Tplot, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path, save_tag="rx" , fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end


@everywhere function self_aligning_step(rx_step,J, Dr,seed, save_folder_path; Tsave=100, Tplot=nothing)

    external_forces = ( ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,J),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = [soft_disk_force([1, 2],[1 1; 1 1])]

    local_dofevolvers = [overdamped_xvf_evolver(1),overdamped_pq_evolver(1)]
    global_dofevolvers = []
    field_dofevolvers = []

    sizes = rx_step.system.sizes
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
        p_i.v0[1] = 0.01

    end
    system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers,field_dofevolvers,false,rx_step.system.rcut_pair_global);

    #Run integration
    sim = Euler_integrator(system,1e-2, 5e3,Tsave=Tsave,seed=seed, Tplot=Tplot, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path, save_tag="sa" , fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim
end

rx_result= relaxation_step("",Tsave=nothing, Tplot=100)
self_aligning_step(rx_result,0.8, 0.01,1, ""; Tsave=nothing, Tplot=100);


Drs = [0., 0.01, 0.1, 1, 10] 
Js=[0.1, 1.]

seeds = reshape( collect(1:length(Drs)*length(Js)), (length(Drs),length(Js)) )

for j in eachindex(Js)
    @sync @distributed for i in eachindex(Drs)

        J = Js[j]
        Dr = Drs[i]

        seed = seeds[i,j]

        display("Running")

        save_folder_path = joinpath(homedir(),"sa","vary_J_Dr_largeN","simdata", "J_$J","Dr_$Dr","seed_$seed");
        print(save_folder_path)

        rx_step = relaxation_step(save_folder_path)

        sim = self_aligning_step(rx_step,J,Dr,seed, save_folder_path);
    end
end