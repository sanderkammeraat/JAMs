include("../src/Engine.jl")
using Random, Distributions

function simulation()

    external_forces = [ABP_2d_angular_noise()]
    karray = @SMatrix [1.0 1.0; 1.0 1.0]

    ϵarray = @SMatrix [0.01 0.01; 0.01 0.1]
    pair_forces = (pairABP_force(1.1), soft_atre_type_force(karray,ϵarray ))
    #pair_forces = (pairABP_force(1.1), soft_disk_force())
    
    
    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N1 = 100
    N2 = 300
    L=100.
    v0 = 0.8
    Dr = 0.01
    poly = 0.2

    fr = 0.2

    initial_state = [ PolarParticle2dNtype(i,1,rand(Uniform(1-poly, 1+poly)),1.,v0,Dr,rand(Uniform(-L/2*fr, L/2*fr),2),[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],[0.],[0.,0.],1,[0.,0.],[0.,0.]) for i=1:N1];
    for i in N1+1:N1+N2
        push!(initial_state,PolarParticle2dNtype(i,2,rand(Uniform(1-poly, 1+poly)),1.,v0,Dr,rand(Uniform(-L/2*fr, L/2*fr),2),[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],[0.],[0.,0.],1,[0.,0.],[0.,0.]))
    end

    size = [L,L];

    system = System(size, initial_state, external_forces, pair_forces , dofevolvers, true);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    states = Euler_integrator(system, 0.1, 100000, 100000, 10, (plot_type_sized_points!, plot_velocity_vectors!));
    0
end



states = simulation()

@time simulation()

@profview simulation()