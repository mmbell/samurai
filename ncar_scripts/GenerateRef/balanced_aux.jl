# *****************************************************************************
# balanced_aux.jl
# 
# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.0.0
#
# This script will house auxilary functions required to compute the balanced
# vortex following Smith (2006; Tellus)
#
# Function list
# ambientprof
# writeout_ncvars
# *****************************************************************************

#==============================================================================
ambientprof

Compute the ambient profiles for exner function and density potential temperature
** Check units!
height = km
rhoa = kg/m^3
temp = K
theta = K
qv = g/kg
==============================================================================#

function ambientprof(height::Real,dbz::Real,rhoa::Real,temp::Real,theta::Real,qv::Real)

    dbz > 50.0 ? dbz = 50.0 : dbz = dbz
    zz = 10.0^(dbz * 0.1)
    rainmass = (zz / 14630.0)^(0.6905)
    icemass =  (zz / 670.0)^(0.5587)
    mixed_dbz = 20.0
    rain_dbz = 30.0
    hlow = 5.0
    melting_zone = 1.0
    hhi = hlow + melting_zone
    if dbz > mixed_dbz && dbz <= rain_dbz
      weightr = (dbz - mixed_dbz) / (rain_dbz - mixed_dbz)
      weights = 1.0 - weightr
      icemass = (rainmass * weightr + icemass * weights) / (weightr + weights)
    elseif dbz > 30
      icemass = rainmass
    end
    precipmass = rainmass * (hhi - height) / melting_zone +
                 icemass * (height - hlow) / melting_zone
    height < hlow ? precipmass = rainmass : nothing
    height > hhi  ? precipmass = icemass  : nothing
    # Compute the exner function and thetarho
    qr = precipmass/rhoa
    exner = temp/theta
    thetarho = theta * (1 + qv / 1000.0) / (1 + qv / 1000.0 + qr / 1000.0)
    return exner,thetarho
end

#==============================================================================
writencvars

This function will be used to write all the variables required by the
thermodynamic retrieval to a NetCDF file

** Note that this is specific to what's required in the thermodynamic retrieval
   but can be generalized for different purposes
==============================================================================#

function writeout_ncvars(filename::AbstractString,coorddict,datadict::OrderedDict)

    # Define the NcVar object
    ncvars = NcVar[]
    # Define the dimensions
    xdim = NcDim("x",coorddict["x"])
    ydim = NcDim("y",coorddict["y"])
    zdim = NcDim("z",coorddict["z"])
    tdim = NcDim("t",[1.0])
    # Define lon and lat variables independently since they're one-dimensional
    lonvar = NetCDF.NcVar("lon",[xdim],t=Float64)
    push!(ncvars,lonvar);
    latvar = NetCDF.NcVar("lat",[ydim],t=Float64)
    push!(ncvars,latvar);
    # Loop over all the variable names in dictionary and push them to ncvars
    for (index,varname) in enumerate(keys(datadict))
        atts = Dict("long_name" => varname, "missing_value" => -999.0, "_FillValue" => -999.0)
        push!(ncvars,NetCDF.NcVar(varname,[xdim,ydim,zdim,tdim],atts=atts,t=Float64))
    end
    # Create the output NetCDF file
    isfile("./" * filename) ? rm("./" * filename) : nothing
    nc = NetCDF.create(filename,ncvars)
    # Place lon and lat in the NetCDF file first
    NetCDF.putvar(nc,"lon",coorddict["lon"])
    NetCDF.putvar(nc,"lat",coorddict["lat"])
    # Loop over remaining variables and push them to NetCDF file
    for (index,varname) in enumerate(keys(datadict))
        NetCDF.putvar(nc,varname,datadict[varname])
    end
    ncclose(filename)
    return nothing
end
