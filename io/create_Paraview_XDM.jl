using HDF5

#DISCLAIMER: this and only this script has been constructed using generative AI


function write_generic_xdmf_2d(h5filename::String, xdmf_filename::String)
    h5 = h5open(h5filename, "r")

    frame_keys = keys(h5["/frames"])
    frame_indices = sort(parse.(Int, collect(frame_keys)))
    first_frame = frame_indices[1]

    all_keys = keys(h5["/frames/$(first_frame)"]) |> collect

    coord_keys = Set(["x", "y"])
    other_keys = setdiff(all_keys, coord_keys)

    # Group potential vector fields by shared base
    grouped = Dict{String, Vector{String}}()
    for key in other_keys
        base = replace(key, r"[xy]$" => "")
        push!(get!(grouped, base, String[]), key)
    end

    n_particles = length(h5["/frames/$(first_frame)/x"])
    close(h5)

    open(xdmf_filename, "w") do io
        println(io, """<?xml version="1.0" ?>""")
        println(io, """<Xdmf Version="3.0">""")
        println(io, """  <Domain>""")
        println(io, """    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">""")

        zero_str = join(fill("0.0", n_particles), " ")

        for (i, frame) in enumerate(frame_indices)
            println(io, """
      <Grid Name="Frame$frame" GridType="Uniform">
        <Time Value="$(i - 1)"/>""")

            # Open file to check types
            h5 = h5open(h5filename, "r")
            frame_group = h5["/frames/$frame"]

            # Add scalar string and numeric metadata as <Information>
            for (key, dset) in pairs(frame_group)
                if key in coord_keys
                    continue
                end

                try
                    if length(size(dset)) == 0  # Scalar
                        val = dset[]
                        val_str = (val isa String) ? val : string(val)
                        println(io, """        <Information Name="$key" Value="$val_str"/>""")
                    end
                catch e
                    @warn "Skipping $key due to read error: $e"
                end
            end

            println(io, """
        <Topology TopologyType="Polyvertex" NumberOfElements="$n_particles"/>
        <Geometry GeometryType="XYZ">
          <DataItem ItemType="Function" Function="JOIN(\$0, \$1, \$2)" Dimensions="$n_particles 3">
            <DataItem Format="HDF" Dimensions="$n_particles">$h5filename:/frames/$frame/x</DataItem>
            <DataItem Format="HDF" Dimensions="$n_particles">$h5filename:/frames/$frame/y</DataItem>
            <DataItem Format="XML" Dimensions="$n_particles">$zero_str</DataItem>
          </DataItem>
        </Geometry>
""")

            # Add vector and scalar array attributes
            for (base, keys) in grouped
                full_keys = sort(keys)
                dset_shapes = [size(frame_group[k]) for k in full_keys if haskey(frame_group, k)]

                if length(full_keys) == 2 &&
                   all(endswith.(full_keys, ["x", "y"])) &&
                   all(s -> s == (n_particles,), dset_shapes)

                    println(io, """
        <Attribute Name="$base" AttributeType="Vector" Center="Node">
          <DataItem ItemType="Function" Function="JOIN(\$0, \$1, \$2)" Dimensions="$n_particles 3">
            <DataItem Format="HDF" Dimensions="$n_particles">$h5filename:/frames/$frame/$(full_keys[1])</DataItem>
            <DataItem Format="HDF" Dimensions="$n_particles">$h5filename:/frames/$frame/$(full_keys[2])</DataItem>
            <DataItem Format="XML" Dimensions="$n_particles">$zero_str</DataItem>
          </DataItem>
        </Attribute>
""")
                else
                    for key in full_keys
                        if !haskey(frame_group, key)
                            continue
                        end
                        dset = frame_group[key]
                        if size(dset) == (n_particles,)
                            println(io, """
        <Attribute Name="$key" AttributeType="Scalar" Center="Node">
          <DataItem Format="HDF" Dimensions="$n_particles">$h5filename:/frames/$frame/$key</DataItem>
        </Attribute>
""")
                        end
                    end
                end
            end

            close(h5)
            println(io, "      </Grid>")
        end

        println(io, """    </Grid>""")
        println(io, """  </Domain>""")
        println(io, """</Xdmf>""")
    end

    println("✅ XDMF written to: $xdmf_filename")
end





# Base folder

base_folder_path = "/Volumes/T7_Shield/sa/survey/hex_disordered/phi_1/Nlin_20/t5e4/simdata/J_0.1/Dr_0.01/seed_1/"

raw_data_file_name = "sa_raw_data.h5"

write_generic_xdmf_2d(joinpath(base_folder_path,raw_data_file_name), joinpath(base_folder_path,"particles.xdmf"))


#
base_folder_path = "/Users/kammeraat/mounting/pi-henkes/kammeraat/pairAN/sim_for_hdf5reader/simdata/"

raw_data_file_name = "raw_data.h5"

write_generic_xdmf_2d(joinpath(base_folder_path,raw_data_file_name), joinpath(base_folder_path,"particles.xdmf"))

