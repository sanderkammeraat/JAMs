

using Random, Distributions
using LinearAlgebra
using StaticArrays


#abstract type Force end

abstract type ExternalForce end

abstract type FieldForce end

abstract type PairForce end

#Neighbour forces: forces where the contributions of other particles needs to be known to take e.g. average 
abstract type NeighbourhoodForce end 
#If particle overlap or special range smaller than rcut_pair_global needs to be taken into account this should be done in the contribute function.
#The engine part will only check against dxn<=rcut_pair_global, so any special condition to strictly sort neighbours needs to in contribute function.

include(joinpath("Forces","ExternalForces.jl"))
include(joinpath("Forces","FieldForces.jl"))
include(joinpath("Forces","PairForces.jl"))
include(joinpath("Forces","NeighbourhoodForces.jl"))
