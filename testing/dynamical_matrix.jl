
include("../analysis/AnalysisPipeline.jl")


raw_data_file = load_file("/Users/kammeraat/seed_1/sa_raw_data.h5")

support_raw_data_file = load_file("/Users/kammeraat/seed_1/ra_raw_data.h5")

frames = raw_data_file["frames"]

system = raw_data_file["system"]

frames_support = support_raw_data_file["frames"]

x0 = frames_support[string(length(frames_support))]["x"]
y0 = frames_support[string(length(frames_support))]["y"]

k = system["forces"]["pair"]["soft_disk_force"]["karray"]

J = system["forces"]["external"]["self_align_with_v_unit_force"]["β"]

R = frames["1"]["R"]

type = frames["1"]["type"]

Dr = frames["1"]["Dr"][1]

begin
D_wb = construct_D(x0, y0, k, R, type)

interior_indmax = sum(type.==1)

D = Symmetric(D_wb[1:2*interior_indmax,1:2*interior_indmax])

eigenmodes = diagonalize_D(D)
end
begin
M=zeros(2*length(x0), 2*length(x0))
Id = [1 0 ; 0 1]
nij = zeros(2,2)

for i in eachindex(x0)
    for j in eachindex(x0)
        if i!==j

            dx = x0[j] - x0[i]
            dy = y0[j] - y0[i]

            dr2 = dx^2 + dy^2
            dr = sqrt(dr2)
            kij = k[type[i],type[j]]

            Rij = R[i]+R[j]
            if dr<Rij
                nij[1,1] =  dx*dx/dr2
                nij[1,2] =  dx*dy/dr2
                nij[2,1] =  dx*dy/dr2
                nij[2,2] =  dy*dy/dr2
                M[2i-1:2i,2j-1:2j] .= -kij .* nij .-kij .* (1 - Rij/dr) .* (Id .- nij)
            end

        end

    end

end
for i in eachindex(x0)
    for j in eachindex(x0)
        if i!==j
            M[2i-1:2i,2i-1:2i].+= -M[2i-1:2i,2j-1:2j]
        end
    end
end
end

eigenmodes_M = diagonalize_D(M[1:2*interior_indmax,1:2*interior_indmax])

eigenmodes_D = diagonalize_D(D[1:2*interior_indmax,1:2*interior_indmax])

eigenmodes_M["eigvecs"][:,2]
eigenmodes_D["eigvecs"][:,2]
#CHECK ✓ it works

close(raw_data_file)
close(support_raw_data_file)


#%%%



using Roots




Nint = interior_indmax



eigvals = eigenmodes_D["eigvals"]

function root_a(x)


    tau = 1/0.1

    J = 0.01


    A = sqrt.((1/tau + J * x) .* eigvals[2:end])




    B = @. eigvals[2:end] - J/x + J*x + 1/tau

    numerator =  @. B - sqrt( B^2 - 4*A^2) 

    denominator =  @. B * sqrt(2*B * (numerator) - 4*A^2 )

    return   @. 2/tau * pi/2 * sum( numerator/denominator  )/2 - Nint * x^2



end

tau = 1/0.1
sqrt(1/Nint*sum(1 ./(2 .+ 2*tau .* eigvals)))


root_a(0.0001)

find_zero(root_a, 0.01)
x=0.9
@.eigvals[2:end] - J/x + J*x + 1/tau