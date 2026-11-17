##---------------------------------------------------------------
## 05_clean_trait_data.jl: Convert TRY and XFT trait data to usable format
##
## author: jacob levine; jacob.levine@utah.edu
##---------------------------------------------------------------

include(joinpath(@__DIR__, "..", "src", "julia", "utils.jl"))

using CSV, DataFrames, DataFramesMeta
using Statistics
using LightGraphs
using DelimitedFiles
using Plots

##---------------------------------------------------------------
## for hydraulic traits, preferentially use XFT values, as these
## can be cleaned for method
##---------------------------------------------------------------

species_list = CSV.read(datadir("fulldata_specieslist.csv"), DataFrame)

xft_traits = CSV.read(datadir("xft", "xylem_traits.csv"), DataFrame)

# Rename columns to avoid issues with "." in column names
rename!(xft_traits, "Cleaned.binomial" => :cleaned_binomial, "Plant.organ" => :plant_organ,
        "P50.method" => :p50_method, "P50..MPa." => :p50_mpa)

## filter out all species not found in the FIA database
xft_traits = filter(row -> row.cleaned_binomial in species_list.species, xft_traits)

# Filter invalid plant organ and P50 methods
xft_traits = filter(row -> !ismissing(row.plant_organ) && !ismissing(row.p50_method) &&
    row.plant_organ in ["S"] &&
    row.p50_method in ["DH", "CE", "CA", "AD", "AS", "Pn"], xft_traits)

# Convert P50 values to Float64 and handle missing values
xft_traits.p50_mpa .= [isnothing(s) ? missing : try parse(Float64, s) catch e missing end for s in xft_traits.p50_mpa];

## remove unrealistic p50 values
xft_traits = filter(row -> !ismissing(row.p50_mpa) && row.p50_mpa < -0.5, xft_traits);

## group by species and calculate mean P50 value
p50_means_xft = combine(groupby(xft_traits, [:cleaned_binomial]), :p50_mpa => mean => :p50_mean);
p50_means_xft = p50_means_xft[completecases(p50_means_xft),:];

##---------------------------------------------------------------
## Now pull remaining p50 data available from TRY
##---------------------------------------------------------------
try_traits = CSV.read(datadir("try", "try_traits.txt"), DataFrame);
tanoak = CSV.read(datadir("try", "tanoak_traits.txt", DataFrame), DataFrame);
try_traits = vcat(try_traits, tanoak);

## Check if file exists, and write it if it doesn't
if !isfile(datadir("try", "trait_traits.csv"))
    CSV.write(datadir("try", "trait_traits.csv"), try_traits)
end

## select only relevant columns
select!(try_traits, Not([:SpeciesName, :AccSpeciesID, :ErrorRisk, :StdValueStr, :Column29]));

## filter for p50 traits
p50_try = filter(row -> !ismissing(row.TraitName) && row.TraitID == 719, try_traits);

## pull out only P50 values, data also includes P88, P12, HSM, etc
acc_names = ["Stem P50", "Xylem water potential at which 50% of conductivity is lost (P50)",
             "Mean P50 including all data"];
p50_try = filter(row -> row.DataName ∈ acc_names, p50_try);

## remove duplicate data
p50_try = p50_try[[!(x in ["Xylem Functional Traits (XFT) Database: Nature Subset",
                           "Xylem Functional Traits (XFT) Database"]) for x in p50_try.Dataset],:];

## remove species already in XFT data
p50_try = p50_try[[!(x in p50_means_xft.cleaned_binomial) for x in p50_try.AccSpeciesName],:];

## remove missing observations
p50_try = p50_try[.!ismissing.(p50_try.StdValue),:];

## correct data entry errors and filter unrealistic values
p50_try[p50_try.StdValue .> 0, :StdValue] .= -1 .* p50_try[p50_try.StdValue .> 0, :StdValue];
filter!(row -> row.StdValue < -0.5, p50_try);

p50_means_try = combine(groupby(p50_try, :AccSpeciesName), :StdValue => (x -> mean(skipmissing(x))) => :P50_mean)


## merge p50 data from both sources and label the source database
rename!(p50_means_try, :AccSpeciesName => :species, :P50_mean => :p50_mean)
rename!(p50_means_xft, :cleaned_binomial => :species)

p50_means_xft.p50_database .= "xft"
p50_means_try.p50_database .= "try"

p50_means = vcat(p50_means_xft, p50_means_try)


##---------------------------------------------------------------
## Amax
##---------------------------------------------------------------

Amax = filter(row -> !ismissing(row.TraitName) && row.TraitID == 53, try_traits);

## filter for area-based Amax only
filter!(row -> row.DataID ∈ [68,2356,1042], Amax);

## remove unrealistic values
filter!(row -> !ismissing(row.StdValue) && row.StdValue .< 40, Amax)

## aggregate Amax data
Amax_means = combine(groupby(Amax, :AccSpeciesName), :StdValue => (x -> mean(skipmissing(x))) => :Amax_mean);
Amax_means.Amax_database .= "try";
rename!(Amax_means, :AccSpeciesName => :species);

##---------------------------------------------------------------
## Stem density
##---------------------------------------------------------------

## get wood density (stem specific density; SSD)
ssd = filter(row -> !ismissing(row.TraitName) &&
    row.TraitID == 4, try_traits);

## filter for relevant metrics
acc_names = ["Wood density; stem specific density; wood specific gravity (SSD)",
             "Wood density; stem specific density; wood specific gravity",
             "Wood density after drying at 100C",
             "Wood density at 0% humidity"];
filter!(row -> row.DataName in acc_names, ssd);

## remove unrealistic values
filter!(row -> row.StdValue < 0.8, ssd)

## aggregate stem density data
ssd_means = combine(groupby(ssd, :AccSpeciesName), :StdValue => (x -> mean(skipmissing(x))) => :ssd_mean);
ssd_means.ssd_database .= "try";
rename!(ssd_means, :AccSpeciesName => :species)

##---------------------------------------------------------------
## Combine data and write to file
##---------------------------------------------------------------

full_traits = outerjoin(p50_means, Amax_means, ssd_means, on = :species);

## reported value of like 250!
full_traits[full_traits.species .== "Ailanthus altissima", :Amax_mean] .= missing;

## sort and write to file
sort!(full_traits, :species);
CSV.write(datadir("full_traits.csv"), full_traits)






