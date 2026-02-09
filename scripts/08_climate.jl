##-------------------------------------------------------------------
## 08_climate_data.jl Pull and process climate data
##
## author: jacob levine; jacob.levine@utah.edu
##-------------------------------------------------------------------

include(joinpath(@__DIR__, "..", "src", "julia", "utils.jl"))

using CSV, DataFrames
using NCDatasets
using ProgressBars, BenchmarkTools
using Plots
using Statistics
using Downloads, Printf, Serialization

# Load 0.25 degree grid centers
data = CSV.read(datadir("templates", "grid_cells.csv"), DataFrame)
latlon = unique(data[:, [:LON, :LAT]])
latlon.grid_cell_id = data.grid_cell_id
latlon = sort(latlon, :LON)
npoints = nrow(latlon)

# Load Terraclimate 4km lon/lat
ds = NCDatasets.Dataset(datadir("terraclim", "TerraClimate_ppt_1958.nc"))
lon = ds["lon"][:]
lat = ds["lat"][:]
close(ds)

##---------------------------------------------------------------
## 00. Precompute 4km pixel overlaps with 0.25° grid
##---------------------------------------------------------------

# Function to get bounding box of 0.25-degree grid
function build_grid_cell_polygon(lon::Float64, lat::Float64, res::Float64=0.25)
    half = res / 2.0
    x = [lon - half, lon + half, lon + half, lon - half]
    y = [lat - half, lat - half, lat + half, lat + half]
    return (xmin=minimum(x), xmax=maximum(x), ymin=minimum(y), ymax=maximum(y))
end

# Precompute pixel indices for each 0.25° grid cell using bounding box filter
grid_to_terracell = Vector{Vector{CartesianIndex{2}}}(undef, npoints)

for i in ProgressBar(1:npoints)
    lon_c, lat_c = latlon.LON[i], latlon.LAT[i]
    bbox = build_grid_cell_polygon(lon_c, lat_c)

    lon_inds = findall(x -> (bbox.xmin <= x <= bbox.xmax), lon)
    lat_inds = findall(y -> (bbox.ymin <= y <= bbox.ymax), lat)

    inds = [CartesianIndex(lon_i, lat_j) for lon_i in lon_inds for lat_j in lat_inds]
    grid_to_terracell[i] = inds
end

serialize(datadir("terraclim", "grid_to_terracell_map.ser"))

grid_to_terracell = deserialize(datadir("terraclim", "grid_to_terracell_map.ser"));


##------------------------------------------------------------------
## 01. Yearly data, for determining drought occurrence and intensity
##------------------------------------------------------------------

# Inputs
years = 1995:2023
n_years = length(years)
vars = ["def", "ppt", "PDSI"]
out_file = datadir("terraclim", "terraclim_processed_grid.csv")

# Prepare metadata
npoints = nrow(latlon)
lon_vals = latlon.LON
lat_vals = latlon.LAT

# Open output CSV file and write header
open(out_file, "w") do io
    println(io, "lon,lat,year,def,ppt,pdsi")  # lowercase headers
end

# Helper function: aggregate 3D array over time (months)
function aggregate_monthly(data_3d, v::String)
    nlon, nlat = size(data_3d, 1), size(data_3d, 2)
    derived = Array{Float64}(undef, nlon, nlat)
    @inbounds for xi in 1:nlon, yi in 1:nlat
        vals = skipmissing(view(data_3d, xi, yi, :))
        derived[xi, yi] = isempty(vals) ? NaN :
                          (v == "ppt" ? sum(vals) : mean(vals))
    end
    return derived
end

# Main processing loop
for y in ProgressBar(years)
    row_data = DataFrame(lon = lon_vals, lat = lat_vals, year = fill(y, npoints))

    for v in vars
        ds = NCDataset(datadir("terraclim", "TerraClimate_$(v)_$(y).nc"))
        data_3d = ds[v][:,:,:]
        close(ds)

        derived = aggregate_monthly(data_3d, v)

        result = Vector{Float64}(undef, npoints)
        Threads.@threads for i in 1:npoints
            inds = grid_to_terracell[i]
            vals = [derived[I] for I in inds if isfinite(derived[I])]
            result[i] = isempty(vals) ? NaN : mean(vals)
        end

        row_data[!, lowercase(v)] = result
    end

    # Append to CSV
    CSV.write(out_file, row_data, append = true)
end

out = CSV.read(out_file, DataFrame)

out.key = string.(round.(out.lat, digits=4)) .* "_" .* string.(round.(out.lon, digits=4))
latlon.key = string.(round.(latlon.LAT, digits=4)) .* "_" .* string.(round.(latlon.LON, digits=4))

leftjoin!(out, latlon, on = :key)

CSV.write(out_file, out)


##---------------------------------------------------------------
## 02. Climatological data, for determining MAT and MAP
##---------------------------------------------------------------

# Download annual climatological data from TerraClimate
local_dir = datadir("terraclim", "local")
isdir(local_dir) || mkdir(local_dir)

vars = ["tmin", "tmax", "ppt"]
years = 1961:2020

for var in vars, yr in years
    fname = @sprintf("TerraClimate_%s_%d.nc", var, yr)
    url = "http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data/$fname"
    dest = joinpath(local_dir, fname)
    if !isfile(dest)
        @printf("Downloading %s\n", fname)
        try
            Downloads.download(url, dest)
        catch e
            @printf("Failed to download %s: %s\n", fname, e)
        end
    end
end

## function to stack data for a given variable
function stack_data(var)
    files = [joinpath(local_dir, @sprintf("TerraClimate_%s_%d.nc", var, yr)) for yr in years]
    ds0 = NCDataset(files[1])
    data0 = ds0[var][:,:,:]
    lon, lat = ds0["lon"][:], ds0["lat"][:]

    mv = ds0[var].attrib["_FillValue"]

    close(ds0)

    nlon, nlat = size(data0, 1), size(data0, 2)
    nt = length(files)
    stack = Array{Union{Missing, Float64}}(undef, nlon, nlat, nt)

    for (i, f) in enumerate(files)
        ds = NCDataset(f)
        arr = ds[var][:,:,:]
        arr = replace(arr, mv => missing)
        if var == "ppt"
            arr_agg = sum(arr, dims=3)  # sum precipitation
        else
            arr_agg = mean(arr, dims=3)  # mean temperature or PDSI
        end
        close(ds)
        stack[:, :, i] = arr_agg
    end

    return lon, lat, stack
end

## stack data for tmin, tmax, and ppt
lon, lat, tmin_stack = stack_data("tmin")
_, _, tmax_stack = stack_data("tmax")
_, _, ppt_stack  = stack_data("ppt")

## Compute means across years
function compute_mean(stack)
    nlon, nlat, _ = size(stack)
    out = Array{Union{Missing, Float64}}(undef, nlon, nlat)
    for i in 1:nlon, j in 1:nlat
        vals = skipmissing(view(stack, i, j, :))
        out[i, j] = isempty(vals) ? missing : Float64(mean(vals))
    end
    return out
end

tmin_mean = compute_mean(tmin_stack)
tmax_mean = compute_mean(tmax_stack)
ppt_mean = compute_mean(ppt_stack)

## Save climatological data
serialize(datadir("terraclim", "tmin_mean.csv"), tmin_mean)
serialize(datadir("terraclim", "tmax_mean.csv"), tmax_mean)
serialize(datadir("terraclim", "ppt_mean.csv"), ppt_mean)

##---------------------------------------------------------------
## Attach to gridcells
##---------------------------------------------------------------

tmin_mean = deserialize(datadir("terraclim", "tmin_mean.ser"))
tmax_mean = deserialize(datadir("terraclim", "tmax_mean.ser"))
ppt_total = deserialize(datadir("terraclim", "ppt_mean.ser"))

npoints = nrow(latlon)

## function to compute mean across indices, handling missing values
function mean_at_indices(arr, inds)
    vals = [arr[I] for I in inds if !ismissing(arr[I]) && !isnan(arr[I])]
    isempty(vals) ? missing : mean(vals)
end

# Initialize result vectors
mat = Vector{Union{Missing, Float64}}(undef, npoints)
map = Vector{Union{Missing, Float64}}(undef, npoints)

for i in ProgressBar(1:npoints)
    indices = grid_to_terracell[i]

    # Extract and average climatologies over the grid cell
    tmin = mean_at_indices(tmin_mean, indices)
    tmax = mean_at_indices(tmax_mean, indices)
    map[i]  = mean_at_indices(ppt_mean, indices)

    # Compute mean annual temperature (MAT) as the mean of tmin and tmax
    mat[i] = mean(skipmissing((tmin, tmax)))

end

df_out = DataFrame(
    grid_cell_id = latlon.grid_cell_id,
    lon = latlon.LON,
    lat = latlon.LAT,
    MAT = mat,
    MAP = map,
)

CSV.write(datadir("terraclim", "climatologies_grid.csv"), df_out)
