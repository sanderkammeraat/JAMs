


include("../analysis/AnalysisPipeline.jl")

#base_folder = "/Users/kammeraat/mounting/alicedata1_kammeraatsc1/sa/statistics/hex_disordered/phi_1.3/Nlin_20"

filepath = "/Users/kammeraat/mounting/data2_kammeraat/sa/statistics/hex_disordered/phi_1.3/Nlin_20/analysis/v0_0.01/Dr_0.1/J_0.01/seed_1.h5"


file = load_file(filepath)


keys(file)


v_projs = file["projs"]["v_projs"]


begin
    
    f = Figure()
    ax = Axis(f[1,1])

    scatterlines!( v_projs[1,:])

    vrms = file["vrms_particle_avg_time"]

    scatterlines!(vrms/0.01)

    display(f)
end

#Conclusion: indeed steady state holds below and above the transition: ✓

close(file)