include("../src/Engine.jl")
using Random, Distributions

function simulation()

    external_forces = (ABP_2d_propulsion_force(), ABP_2d_angular_noise())

    karray = @SMatrix [1.0 1.0; 1.0 1.0]

    ϵarray = @SMatrix [0.5 0.5; 0.5 0.5]

    pair_forces = [soft_atre_type_force(karray,ϵarray )]

    
    
    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N1=500
    N2 = 20
    ϕ = 0.6
    L=sqrt((N1+N2)*pi/ϕ)
    poly = 0.0000002

    Dr = 0.01
    v01 =0.05
    v02 = 0.2
    fr=0.2


    initial_state = [ PolarParticle2dNtype(i,1,rand(Uniform(1-poly, 1+poly)),1.,v01,Dr,rand(Uniform(-L/2*fr, L/2*fr),2),[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],[0.],[0.,0.],1,[0.,0.],[0.,0.]) for i=1:N1];
    for i in N1+1:N1+N2
        push!(initial_state,PolarParticle2dNtype(i,2,rand(Uniform(1-poly, 1+poly)),1.,v02,Dr,rand(Uniform(-L/2*fr, L/2*fr),2),[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],[0.],[0.,0.],1,[0.,0.],[0.,0.]))
    end

    size = [L,L];

    system = System(size, initial_state, external_forces, pair_forces , dofevolvers, false);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    states = Euler_integrator(system, 0.1, 100000, 100000, 10, [plot_type_sized_points!, plot_directors!]);
    return states

end


states = simulation()

