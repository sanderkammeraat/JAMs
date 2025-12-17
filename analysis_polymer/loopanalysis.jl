include("analysisfunctions.jl")
include("loaddata.jl")

ksd = [1.]
kbend = [.3]
kstretch = [1.]
fstretch = [.7]
p = [0.01, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2]
kperp = [0]
kpar = [-1]
Npol = [150]
N = [1500]

analysis = [MSD, average_velocity]

for ksd_value in ksd
for kbend_value in kbend
for kstretch_value in kstretch
for fstretch_value in f_stretch
for p_value in p
for kperp_value in kperp
for kpar_value in kpar
for Npol_value in Npol
for N_value in N 
    for function in analysis
        function()
    end
end
end
end
end
end
end
end
end
end

