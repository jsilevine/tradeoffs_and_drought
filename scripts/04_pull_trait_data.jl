##---------------------------------------------------------------
## 04_pull_trait_data.jl: Generate TRY species code list for pulling trait data
##
## Author: Jacob Levine; jacob.levine@utah.edu
##---------------------------------------------------------------

include(joinpath(@__DIR__, "..", "src", "julia", "utils.jl"))

using CSV, DataFrames, DataFramesMeta
using Statistics, LightGraphs
using DelimitedFiles
using Plots

## load data
states = CSV.read(datadir("state_list.csv"), DataFrame)
try_list = CSV.read(datadir("try_specieslist.txt"), DataFrame)
fia_list = CSV.read(datadir("fia_specieslist.csv"), DataFrame)

## clean up fia_list to match FIA tree data
select!(fia_list, ["FIA Code", "Common Name", "Genus", "Species", "PLANTS Code"])
rename!(fia_list, "FIA Code" => :SPCD, "Common Name" => :COMMON_NAME, "Genus" => :GENUS,
       "Species" => :SPECIES, "PLANTS Code" => :PLANTS_CODE)

## create new column with complete species names to match with TRY data
fia_list.species .= string.(fia_list.GENUS) .* " " .* string.(fia_list.SPECIES)

## rename tanoak, bc USFS is living in the 90s
fia_list[fia_list.species .== "Lithocarpus densiflorus", :species] .= "Notholithocarpus densiflorus"

## create empty vector of TRY IDs, to populate from state-level tree data
try_ids = DataFrame(try_id = Vector{Int64}(undef, 0),
                    species = Vector{String}(undef, 0),
                    state = Vector{String}(undef, 0))

## loop through states and pull the TRY species codes that match each unique species
for state in states.state

    println("starting: " * state)

    ## read in state-level tree data
    tree_data = CSV.read(datadir("fia", state, state * "_TREE.csv"), DataFrame)

    ## convert species code to integer for better matching
    tree_data.SPCD .= Int.(tree_data.SPCD)

    ## join FIA species data to state-level tree data
    tree_data = leftjoin(tree_data, fia_list, on = :SPCD)

    ## pull out all unique species from state
    species_list = unique(tree_data.species)

    species_list = species_list[.!ismissing.(species_list), :]

    ## remove entries that don't match FIA codes
    species_list = species_list[.!occursin.("missing", species_list) .&& .!occursin.("unknown", species_list)]

    ## get TRY IDs
    ndata = DataFrame(try_id = filter(row -> row.AccSpeciesName ∈ species_list, try_list).AccSpeciesID,
                      species = filter(row -> row.AccSpeciesName ∈ species_list, try_list).AccSpeciesName,
                      state = repeat([state], length(filter(row -> row.AccSpeciesName ∈ species_list, try_list).AccSpeciesID)))

    try_ids = vcat(try_ids, ndata)

end

try_ids = try_ids[.!nonunique(try_ids[:,[:try_id, :species]]),[:try_id, :species]]
fia_list = fia_list[.!nonunique(fia_list[:,[:SPCD]]),:]

try_ids = leftjoin(try_ids, fia_list[:,[:species, :SPCD]], on = :species)

rename!(try_ids, :SPCD => :fia_id)

## write species list to csv for reference
CSV.write(datadir("fulldata_specieslist.csv"), try_ids)

## create .txt file with list of comma seperated TRY species IDs.
## Will paste this directly into TRY data request
open(datadir("try", "species_topull.txt"), "w") do f
    finalid = unique(try_ids.try_id)[length(unique(try_ids.try_id))]
    for id in unique(try_ids.try_id)
        if id == finalid
            print(f, id)
        else
            print(f, string(id) * ",")
        end
        flush(f)
    end
end

## trait list:
## stem specific density / wood density : 4
## xylem hydraulic vulnerability curve: 3479
## xylem hydraulic vulnerability: 719
## photosynthesis rate per leaf area: 53
## A/Ci curve: photosynthetic rate per leaf area: 3364
## Plant relative growth rate: 77
## Plant growth rate: 587


