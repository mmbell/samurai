# *****************************************************************************
# balanced_vortex.jl
#
# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.0.0
#
# This script will compute the balanced state of a vortex following
# Smith (2006; Tellus).
#
# Description of variables to be modified by user
#
# fin, fout
# Input file name (fin), output filename (fout) if user specifes
# writeout = true (see below)
#
# mode
# "hydrostatic" - create a hydrostatic (horizontally homogeneous) reference
#                 state based on ambient profile
# "vortex_cart" - compute the balanced state (gradient wind and hydrostatic)
#                 of a vortex initially on a Cartesian coordinate grid
# "vortex_cyln" - compute the balanced state (gradient wind and hydrostatic)
#                 of a vortex initially on a cylindrical coordinate grid
#
# cenfile, cx, cy
# Specify the center of the storm for a gridded Cartesian analysis
# Only applies to mode = "vortex_cart"
# If center is at (0,0) manually override and make cx = 0, cy = 0
#
# paramBL
# Parameterize the boundary layer when computing the balanced vortex?
# Only applies to mode = "vortex_cart" and mode = "vortex_cylnd"
# Note: Current parameterization simply sets winds below 2 km altitude to
#       the value at 2 km altitude
#
# fcor
# Specify the Coriolis parameter for the analysis
#
# writeout
# Option to write output to NetCDF file in the format required for running
# the thermodynamic retrieval
# *****************************************************************************

using Interpolations, NetCDF, DataStructures, Statistics, NCDatasets
push!(LOAD_PATH, "/Users/jcdehart/julia/JuliaMet/src") # need to modify for other users
using JuliaMet
include("./balanced_aux.jl")

#=========================== BEGIN MODIFY ====================================#
# Specify the input file (currently, script formulated specific to SAMURAI)
    fin = "./sam/output/samurai_XYZ_analysis.nc"
# Choose which mode to run
    mode = "hydrostatic"
    paramBL = true
# Specify the center location
    #cenfile = ".//center_output.nc"
    #cx = mean(ncread(cenfile,"final_xc"))
    #cy = mean(ncread(cenfile,"final_yc"))
    cx = 20
    cy = 4
# Specify NetCDF output file name in format required by thermodynamic retrieval
    fout = "./beltrami_test.nc"
# Specify the Coriolis parameter for case at hand
    fcor = 0.0 # Coriolis
# Write output to NetCDF file?
    writeout = true
#=========================== END MODIFY ======================================#

#==============================================================================
Read in all required variables for computing the balanced state
==============================================================================#

    balvarnames = ["longitude","latitude","x","y","altitude","DBZ","U","V","T","THETA","QV","RHOA"]
    gradvarnames = ["W","DUDX","DVDX","DWDX","DUDY","DVDY","DWDY","DUDZ","DVDZ","DWDZ"]

    fileIN = Dataset(fin,"r")

    lat = fileIN["latitude"][:,:][:,:,:,1]
    lon = fileIN["longitude"][:,:][:,:,:,1]
    x = fileIN["x"][:,:]
    y = fileIN["y"][:,:]
    z = fileIN["altitude"][:,:]
    dbz = fileIN["DBZ"][:,:][:,:,:,1]
    u = fileIN["U"][:,:][:,:,:,1]
    v = fileIN["V"][:,:][:,:,:,1]
    temp = fileIN["T"][:,:][:,:,:,1]
    theta = fileIN["THETA"][:,:][:,:,:,1]
    qv = fileIN["QV"][:,:][:,:,:,1]
    rhoa = fileIN["RHOA"][:,:][:,:,:,1]

    w = fileIN["W"][:,:][:,:,:,1]
    dudx = fileIN["DUDX"][:,:][:,:,:,1]
    dvdx = fileIN["DVDX"][:,:][:,:,:,1]
    dwdx = fileIN["DWDX"][:,:][:,:,:,1]
    dudy = fileIN["DUDY"][:,:][:,:,:,1]
    dvdy = fileIN["DVDY"][:,:][:,:,:,1]
    dwdy = fileIN["DWDY"][:,:][:,:,:,1]
    dudz = fileIN["DUDZ"][:,:][:,:,:,1]
    dvdz = fileIN["DVDZ"][:,:][:,:,:,1]
    dwdz = fileIN["DWDZ"][:,:][:,:,:,1]

    lat = replace(lat,missing=>NaN)
    lon = replace(lon,missing=>NaN)
    x = replace(x,missing=>NaN)
    y = replace(y,missing=>NaN)
    z = replace(z,missing=>NaN)
    dbz = replace(dbz,missing=>NaN)
    u = replace(u,missing=>NaN)
    v = replace(v,missing=>NaN)
    temp = replace(temp,missing=>NaN)
    theta = replace(theta,missing=>NaN)
    qv = replace(qv,missing=>NaN)
    rhoa = replace(rhoa,missing=>NaN)

    w = replace(w,missing=>NaN)
    dudx = replace(dudx,missing=>NaN)
    dvdx = replace(dvdx,missing=>NaN)
    dwdx = replace(dwdx,missing=>NaN)
    dudy = replace(dudy,missing=>NaN)
    dvdy = replace(dvdy,missing=>NaN)
    dwdy = replace(dwdy,missing=>NaN)
    dudz = replace(dudz,missing=>NaN)
    dvdz = replace(dvdz,missing=>NaN)
    dwdz = replace(dwdz,missing=>NaN)

    #lon,lat,x,y,z,dbz,u,v,temp,theta,qv,rhoa = read_ncvars(fin,balvarnames)
    #w,dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz = read_ncvars(fin,gradvarnames,false)

    # Specify an offset based on length of x-array (avoid boundaries)
    offset = convert(Int64,floor(0.1*length(x)))

    xs = offset + 1
    xe = length(x) - offset
    ys = offset + 1
    ye = length(y) - offset

    xinds = xs:xe
    yinds = ys:ye

    lon = lon[xinds]
    lat = lat[yinds]
    x = x[xinds]
    y = y[yinds]
    dbz = dbz[xinds,yinds,:]
    u = u[xinds,yinds,:]
    v = v[xinds,yinds,:]
    w = w[xinds,yinds,:]
    temp = temp[xinds,yinds,:]
    theta = theta[xinds,yinds,:]
    qv = qv[xinds,yinds,:]
    rhoa = rhoa[xinds,yinds,:]
    dudx = dudx[xinds,yinds,:]
    dvdx = dvdx[xinds,yinds,:]
    dwdx = dwdx[xinds,yinds,:]
    dudy = dudy[xinds,yinds,:]
    dvdy = dvdy[xinds,yinds,:]
    dwdy = dwdy[xinds,yinds,:]
    dudz = dudz[xinds,yinds,:]
    dvdz = dvdz[xinds,yinds,:]
    dwdz = dwdz[xinds,yinds,:]

    # Define constants
    g = 9.81 # Gravity constant
    Cp = 1005. # Specific heat of dry air at constant pressure

    r,phi = xy2rp(cx,cy,x,y)

    if mode == "vortex_cart"

       # Now regrid the data to cylindrical coordinates
       dbz_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,dbz)
       dbzlin_rpz = 10.0.^(dbz_rpz ./ 10)
       u_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,u)
       v_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,v)
       temp_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,temp)
       theta_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,theta)
       qv_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,qv)
       rhoa_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,rhoa)

       # Convert u,v to ur,vt
       ur_rpz,vt_rpz = uv2urvt(phi,u_rpz,v_rpz)

       # Compute the azimuthal means
       azmean_dbzlin = nanmean(dbzlin_rpz,2,true)
       azmean_dbz = 10.0 .* log10.(azmean_dbzlin)
       azmean_temp = nanmean(temp_rpz,2,true)
       azmean_theta = nanmean(theta_rpz,2,true)
       azmean_qv = nanmean(qv_rpz,2,true)
       azmean_rhoa = nanmean(rhoa_rpz,2,true)
       azmean_vt = nanmean(vt_rpz,2,true)

       # Specify the radius to compute the ambient profile
       # Following Annette's convention
       ind_prof = findfirst(isequal(abs(x[1])),r)

       # Compute the ambient profile
       amb_exner = Array{Float64}(undef,length(z))
       amb_thetarho = Array{Float64}(undef,length(z))

       for k in eachindex(z)
           amb_exner[k],amb_thetarho[k] = ambientprof(z[k],azmean_dbz[ind_prof,k],azmean_rhoa[ind_prof,k],azmean_temp[ind_prof,k],azmean_theta[ind_prof,k],azmean_qv[ind_prof,k])
       end

    else
	
       # Compute the azimuthal means
       dbzlin = 10.0.^(dbz ./ 10)
       dommean_dbzlin = nanmean(reshape(dbzlin,:,length(z)),1,true)
       dommean_dbz = 10.0 .* log10.(dommean_dbzlin)
       dommean_temp = nanmean(reshape(temp,:,length(z)),1,true)
       dommean_theta = nanmean(reshape(theta,:,length(z)),1,true)
       dommean_qv = nanmean(reshape(qv,:,length(z)),1,true)
       dommean_rhoa = nanmean(reshape(rhoa,:,length(z)),1,true)
       print(dommean_temp)

       # Compute the ambient profile
       amb_exner = Array{Float64}(undef,length(z))
       amb_thetarho = Array{Float64}(undef,length(z))

       for k in eachindex(z)
           amb_exner[k],amb_thetarho[k] = ambientprof(z[k],dommean_dbz[k],dommean_rhoa[k],dommean_temp[k],dommean_theta[k],dommean_qv[k])
       end

    end
    

    # Convert to SI units prior to integration
    rmet = r * 1e3
    zmet = z * 1e3

    # Define drmet and dzmet
    # FLAG - Needs to be updated to account for stretched grids
    drmet = rmet[2] - rmet[1]
    dzmet = zmet[2] - zmet[1]

    # Retrieve finite difference weights (second-order accurate for stretched grids)
    rwgts = fd_weights(rmet)
    zwgts = fd_weights(zmet)

    # Integrate the hydrostatic equation for the ambient profile
    # This will give the hydrostatic pressure at each altitude
    amb_hexner = Array{Float64}(undef,length(z))

    for k in eachindex(z)
        if k == 1
            amb_hexner[k] = amb_exner[1]
        else
            if any(isnan.(amb_thetarho[1:k]))
                amb_hexner[k] = NaN
            else
                amb_hexner[k] = amb_exner[1] - (g/Cp) * newtoncotes(zmet[1:k], 1 ./ amb_thetarho[1:k])
            end
        end
    end

#==============================================================================
Branch off to the specified mode
==============================================================================#

    # Define arrays for pib, trb, dpibdx, dpibdy

    if mode == "hydrostatic"
       pib = Array{Float64}(undef,size(u))
       trb = Array{Float64}(undef,size(u))
    else
       pib = Array{Float64}(undef,size(azmean_vt))
       trb = Array{Float64}(undef,size(azmean_vt))
    end
    
    pib_xy = similar(u)
    trb_xy = similar(u)

    fill!(pib,NaN)
    fill!(trb,NaN)

#============================== HYDROSTATIC ==================================#

    if mode == "hydrostatic"
        # Use the ambient profiles to create hydrostatic base state
        # Can skip ahead and broadcast to pib_xy, trb_xy
        for j in eachindex(x)
            for i in eachindex(y)
                pib_xy[i,j,:] .= amb_hexner
                trb_xy[i,j,:] .= amb_thetarho
            end
        end
        # Compute dpibdx and dpibdy
        # ** NOTE **
        # dpibdx and dpibdy are required in units of 1/km and conversion to
        # SI units takes place in NetCDF_XYZ.cpp of thermodynamic retrieval
        dpibdx = finite_dx(x,pib_xy)
        dpibdy = finite_dy(y,pib_xy)
        # Replace all NaNs with -999.0 required by thermodynamic retrieval
        pib_xy[isnan.(pib_xy)] .= -999.0
        trb_xy[isnan.(trb_xy)] .= -999.0
        dpibdx[isnan.(dpibdx)] .= -999.0
        dpibdy[isnan.(dpibdy)] .= -999.0
    end

#============================== VORTEX_CART ==================================#

    if mode == "vortex_cart"
        if paramBL
            # Set vt below 2 km to the value at 2 km
            # This forces the boundary layer parameterization
            BLtop = closest_ind(z,2.0)
            azmean_vt[:,1:BLtop] .= azmean_vt[:,BLtop]
        end

        # Create interpolation objects for ambient exner and thetarho profiles
        amb_exner_itp = extrapolate(interpolate((zmet,),amb_exner,Gridded(Linear())),NaN)
        amb_thetarho_itp = extrapolate(interpolate((zmet,),amb_thetarho,Gridded(Linear())),NaN)

        # Compute required variables
        C = azmean_vt.^2 ./ rmet + fcor * azmean_vt
        C[1,:] .= NaN # Undefined at r = 0
        dCdz = finite_dsz(zwgts,C)

        # Compute the balanced pressure and density potential temperature
        for k in 1:length(z)-1
            for i in 1:length(r[1:ind_prof])-1
                z_r = 0.
                lnthetarho = 0.
                nan_flag = false
                for ii in i:length(r[1:ind_prof])
                    # Define vertical interpolation weights along the pressure surface
                    zi = (zmet[k] - zmet[1] + z_r)/dzmet
                    z1 = floor(zi)
                    zw2 = zi - z1
                    zw1 = 1.0 - zw2
                    kk = Int(z1) + 1
                    # Don't allow NaNs along integration path 
                    if isnan(C[ii,kk]) || isnan(C[ii,kk+1]) || isnan(dCdz[ii,kk]) || isnan(dCdz[ii,kk+1])
                        nan_flag = true
                        break
                    else
                        # Integrate along the pressure surface
                        C_r = C[ii,kk] * zw1 + C[ii,kk+1] * zw2
                        dCdz_r = dCdz[ii,kk] * zw1 + dCdz[ii,kk+1] * zw2
                        if ii == i || ii == length(r[1:ind_prof])
                            z_r += 0.5 * C_r/g * drmet
                            lnthetarho += 0.5 * dCdz_r/g * drmet
                        else
                            z_r += C_r/g * drmet
                            lnthetarho += dCdz_r/g * drmet
                        end
                    end
                end
                if nan_flag 
                    pib[i,k] = NaN
                    trb[i,k] = NaN
                else
                    # Define vertical interpolation weights at location of ambient profile
                    zi = (zmet[k] - zmet[1] + z_r)/dzmet 
                    z1 = floor(zi)
                    zw2 = zi - z1
                    zw1 = 1.0 - zw2
                    kk = Int(z1) + 1
                    # Define the balanced exner function and density potential temperature 
                    pib[i,k] = amb_exner[kk] * zw1 + amb_exner[kk+1] * zw2
                    trb[i,k] = exp(log(amb_thetarho[kk] * zw1 + amb_thetarho[kk+1] * zw2) - lnthetarho)
                end
            end # radius
        end # altitude

        # Interpolate pib and trb to Cartesian grid
        for k in eachindex(z)
            # Interpolate to the x-y grid
            pib_xy[:,:,k] .= regrid_pol2cart(cx,cy,r[1:ind_prof],x,y,pib[1:ind_prof,k])
            trb_xy[:,:,k] .= regrid_pol2cart(cx,cy,r[1:ind_prof],x,y,trb[1:ind_prof,k])
        end

        # Compute dpibdx and dpibdy
        # ** NOTE **
        # dpibdx and dpibdy are required in units of 1/km and conversion to
        # SI units takes place in NetCDF_XYZ.cpp of thermodynamic retrieval
        dpibdx = finite_dx(x,pib_xy)
        dpibdy = finite_dy(y,pib_xy)

        # Replace all NaNs with -999.0 required by thermodynamic retrieval
        pib_xy[isnan.(pib_xy)] .= -999.0
        trb_xy[isnan.(trb_xy)] .= -999.0
        dpibdx[isnan.(dpibdx)] .= -999.0
        dpibdy[isnan.(dpibdy)] .= -999.0
    end

#==============================================================================
# Write all the variables to the NetCDF file
# Rename dictionary entries according to what's required by thermodynamic
# retrieval and then push additional variables
==============================================================================#

if writeout

# Set u,v missing values to -999.0

    u[isnan.(u)] .= -999.0
    v[isnan.(v)] .= -999.0

# Place the required coordinate arrays in a dict and variables in another dict
# Need to add the singleton (time) dimension back to variables

    coord_dict = OrderedDict()
    thermo_dict = OrderedDict()

    coord_dict["lon"] = lon
    coord_dict["lat"] = lat
    coord_dict["x"] = x
    coord_dict["y"] = y
    coord_dict["z"] = z

    thermo_dict["u"] = reshape(u,(length(x),length(y),length(z),1))
    thermo_dict["v"] = reshape(v,(length(x),length(y),length(z),1))
    thermo_dict["w"] = reshape(w,(length(x),length(y),length(z),1))
    thermo_dict["dudx"] = reshape(dudx,(length(x),length(y),length(z),1))
    thermo_dict["dvdx"] = reshape(dvdx,(length(x),length(y),length(z),1))
    thermo_dict["dwdx"] = reshape(dwdx,(length(x),length(y),length(z),1))
    thermo_dict["dudy"] = reshape(dudy,(length(x),length(y),length(z),1))
    thermo_dict["dvdy"] = reshape(dvdy,(length(x),length(y),length(z),1))
    thermo_dict["dwdy"] = reshape(dwdy,(length(x),length(y),length(z),1))
    thermo_dict["dudz"] = reshape(dudz,(length(x),length(y),length(z),1))
    thermo_dict["dvdz"] = reshape(dvdz,(length(x),length(y),length(z),1))
    thermo_dict["dwdz"] = reshape(dwdz,(length(x),length(y),length(z),1))
    thermo_dict["trb"] = reshape(trb_xy,(length(x),length(y),length(z),1))
    thermo_dict["pib"] = reshape(pib_xy,(length(x),length(y),length(z),1))
    thermo_dict["dpibdx"] = reshape(dpibdx,(length(x),length(y),length(z),1))
    thermo_dict["dpibdy"] = reshape(dpibdy,(length(x),length(y),length(z),1))

# Define an array of strings for the vars needed by thermodynamic retrieval
# (excluding lon,lat,x,y,z coordinate arrays)

    varnames = ["u","v","w","dudx","dvdx","dwdx","dudy","dvdy","dwdy",
                "dudz","dvdz","dwdz","trb","pib","dpibdx","dpibdy"]

# Write to output NetCDF file

    writeout_ncvars(fout,coord_dict,thermo_dict)

end
