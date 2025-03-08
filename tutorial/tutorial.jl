# In this script, we will set up a simple simulation with ABPs and a soft repulsive interaction.
# We will make a live plot during the simulation while saving the simulation data to a folder. 
# After the simulatin is done, we read in this data to show how to acces it.

# Note: as of now, the code is not structured as a package with modules. This means that we will simply include the necessary files to make the custom structs and functions available.

# Let's begin. We load in the the /src/Engine.jl file wich will subsequently load in the other files with particle definitions etc.:
include(joinpath("..","src","Engine.jl"))

# The basic idea of the JAMs workflow is to setup a (physical) system  and pass it to the integrator.
# In JAMs the system should be a instance of the struct System (defined in Engine.jl).

# Let's find out what System needs as input. In the REPL terminal type: ?System or hover with the mouse cursor over System in the following line
System

# We see that we need to set the linear sizes of the system. Here we will first set the number of particles and a packing fraction which will determine the system size

N=100 #number of particles
ϕ = 0.8 #packing fraction
r =1 #radius of particles
L=sqrt(N*pi*r^2/ϕ)

# We do a 2d system, so we define the sizes as 
sizes = [L, L]

# We will not include any fields and updaters related to them so we simply set an empty array for them.
initial_field_state=[]
field_forces = []
field_updaters = []

# We now prepare the initial particle state. This should an array containing structs.
# For the 2d ABPs we will be using the PolarParticle2d struct.
# We can use Julia's list apprehension to fill the initial state of N particles.
# Do in REPL: ?PolarParticle2d to find the meaning of the fields (or hover with mouse above PolarParticle2d in VS Code to get the same information)
# Using rand(Uniform(-L/2,L/2)) we get a random number drawn from a uniform distribution to place the particle randomly in the system
# Using rand(Uniform(-L/2,L/2), 2) we get two random numbers in a vector of size 2
# Using rand(Uniform(-pi, pi)) we get a random angle for the director. Note that we put this between [] to allow it to be changed by the dof evolver
# (fields to a struct are inmmutable, but if the field is an array, the elements in the array can be changed)

# Note that the unwrapped coordinates in xuw are set automatically, so we can set them to 0. JAMs automatically sets the initial unwrapped coordinates equal to the initial wrapped coordinates
# The type in a particle struct allows forces to act on a subset of particles.
# !! For as single particle type, type=1 must be used. If multiple types, start from 1 and use the next integer for a new type (e.g. 1 2 3 etc. and not 1 3 4) !!

Dr = 0.01
v0 = 0.3
initial_particle_state = [ PolarParticle2d(i,1,1,v0,Dr,rand(Uniform(-L/2,L/2), 2),[0,0],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],r,1.,[0.,0.],[0.,0.],[0,0]) for i=1:N]

#Now set up the particle forces 
# The forces take in a type (Int64 or array of Int64 if on multiple particle types) determining on what type of particles the force should act on.
# Note that we categorize noise as a force.
external_forces = [ABP_2d_propulsion_force(1), ABP_2d_angular_noise(1)]

#The soft disk force works on all particles of type 1 and has a stiffness of 1.
pair_forces = [soft_disk_force(1,1)]


#We want overdamped update of the degrees of freedom:
dofevolvers =  [overdamped_evolver!]


# Let's make the system periodic (note the uncapitalized first letter on the bool in Julia)
Periodic= true

# And set the cell lists cut off radius to be 2.5 times the radius of the particles
rcut_pair_global = 2.5*r

#What if your pair forces should not be cut off? Set the cut off radius to the largest linear system size.

# We use all this to setup an instance of System
system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, Periodic,rcut_pair_global);

# Now it's time to integrate. For this we use the Euler_integrator, currently the only integrator in JAMs. (Using it with noise technically makes it a Euler-Maruyama integrator)
# Use ?Euler_integrator or hover above it in VS code for information.
Euler_integrator

# The Euler integrator takes in 3 mandatory arguments: the system as defined above, the Euler time step dt, and the end time t_stop.
# Without any addtional arguments, the initial state in the system will be integrated and the Euler_integrator returns a struct called SIM, which has as fields the final particle and field states, the system, dt, and t_stop.

sim = Euler_integrator(system, 0.1, 1)
#To inspect the final particle_state, do 
sim.final_particle_state


#Typically, we  actuallyt want to save data, so we set the optional argument Tsave: being the save interval in terms of integrator steps (t is advanced by dt every step).
Tsave = 10 # Save every 10 dt steps

# We supply a list of save_functions to collect the data we want to store each time step.
save_functions=[save_2d_polar_θ!]

#Also, we need to set a save_folder_path. It should end with a slash and  if the directory does not exist, it will automatically be created
save_folder_path = "/Users/kammeraat/JAMs_tutorial/"

# If we also want to plot in realtime the system, we need to supply even more optional arguments:

Tplot = 10 # plot every 10 dt time steps
fps = 60 #wait 1/60 before plotting the next frame
plot_functions = [plot_disks!, plot_directors!, plot_velocity_vectors!] # plotting functions to be used

#Let's integrate again, but now with saving and live plotting
sim = Euler_integrator(system, 0.1, 100, Tsave=Tsave, save_functions=save_functions, save_folder_path=save_folder_path, Tplot=Tplot, plot_functions=plot_functions, fps = fps)

############### LOADING data
# The save folder should now contain two files. One is a JAMs_container. Loading it in Julia returns the Julia objects.
# Let's load it (read-only mode "r")
JAMs_container = jldopen(save_folder_path*"JAMs_container.jld2","r")

#We can get the system struct from it (note that system includes the initial state)  by using:
loaded_system = JAMs_container["system"]

#and also the final state
final_particle_state = JAMs_container["final_particle_state"]
#We can use these to setup a sequel simulation!

#The data collected by the save_functions is stored in another file: raw_data.jld2. Which is also readable in e.g. python as it adheres to the hdf5 standard.

raw_data = jldopen(save_folder_path*"raw_data.jld2","r")

# As an example, let's get the x-coordinates of the particles in stored frame number 20. These coordinates are wrapped. 
x_20 = raw_data["frames"]["20"]["x"]

# If we want the unwrapped x-coordinates of the particles, we do
x_20_uw = raw_data["frames"]["20"]["xuw"]



# What is the actual time of this frame? (It is not 20! This is because frame 1 contains the info of t=0).
t_20 = raw_data["frames"]["20"]["t"]

#Note that the raw_data file also contains metadata, e.g. the forces and the field values used:
soft_disk_stiffness = raw_data["system"]["forces"]["pair_forces"]["soft_disk_force"]["karray"]

#Note that by current design, Dr and v0 are saved per particle
v0s = raw_data["frames"]["1"]["v0"]





















s