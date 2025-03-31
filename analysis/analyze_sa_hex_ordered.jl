

include(joinpath("..","src","Engine.jl"))
include("AnalysisFunctions.jl")

raw_data_base_folder = joinpath(homedir(),"sa", "varyJ")

#Make tree to navigate simulation data files structure
tree = Dict()

#Base folders (first parameter axis)
dirs = readdir_filt(raw_data_base_folder)

#Subfolder (second parameter axis)
for dir in dirs
    tree[dir]=Dict()
    subdirs = readdir_filt(joinpath( raw_data_base_folder, dir ))

    for subdir in  subdirs
        tree[dir][subdir]=Dict()

        seeddirs =  readdir_filt(joinpath( raw_data_base_folder, dir, subdir ))

        for seeddir in seeddirs
            tree[dir][subdir][seeddir]= joinpath( raw_data_base_folder, dir, subdir, seeddir )
        end
    end
    

end



