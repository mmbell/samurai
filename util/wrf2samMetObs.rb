=begin
Created on Thr Apr 18 17:47:18 2013
@author: Annette Foerster (foerster@hawaii.edu)
Modified by Ting-Yu Cha (tingyu@colostate.edu) on Fri Feb 4 16:38 2022
=end

require 'numru/netcdf'
require 'date'
require 'geokit'

include NumRu

Geokit::default_units = :kms

def calc_dist(lat1,lon1,lat2,lon2)

  position1 = Geokit::LatLng.new(lat1,lon1)
  position2 = Geokit::LatLng.new(lat2,lon2)

  distance = position1.distance_to(position2)
  return distance
end

# ------------------  TO MODIFY ------------------
center_lat = 40.44625
center_lon = -104.63708
maxradius = 50  # in km
dx = 3*1000 # in m
dy = 3*1000 # in m

init_year = 2021; init_mon  = 7; init_day  = 30; init_hr   = 12
year = 2021; mon  = 7; day  = 30; hr   = 12

filein = "/Users/tingyu/Dropbox/CSU/Research/PhD/samurai/wrf_enkf_output_d02_001.nc"
fileout = "../../wrf_20210730/data/Observation.txt"
fileterrain = "../../wrf_20210730/data/terrain.hgt"

# ------------------  MODIFY END ------------------

# Convert time to seconds and record as time
init_time = DateTime.new(init_year,init_mon,init_day,init_hr)
time = DateTime.new(year,mon,day,hr)
unixtime = time.to_time.to_i.to_s
timestep = ((time-init_time)*24).to_i
time_input = init_time.strftime('%Y-%m-%d_%H:%M:%S')

# Read nc-file
puts "Reading NetCDF for " + init_time.strftime('%Y-%m-%d_%H:%M:%S') + " at timestep " + timestep.to_s + " ..."
puts "Reading NetCDF for " + time.to_s + " at timestep " + timestep.to_s + " ..."
file = NetCDF.open(filein)

raw_data = Hash.new()
adjusted_data = Hash.new()
terrain_data = Hash.new()
dims = ["XTIME","XLAT","XLONG"]
vars = ["U","V","W","T","PH","PHB","P","PB","QVAPOR","REFL_10CM"]
sfcvars = ["U10","V10","T2","PSFC","Q2","HGT"]

dims.each do |dim|
  raw_data[dim] = file.var(dim).get
end

vars.each do |var|
  raw_data[var] = file.var(var).get[true,true,true,timestep]
end

sfcvars.each do |var|
  raw_data[var] = file.var(var).get
end

# Calculate height, total pressure, dry air density and temperature
# Change units of q_vapor from kg/kg to g/kg
puts "Calculating additional variables ..."

size = raw_data["T"].size
rank = raw_data["T"].rank
shape = raw_data["T"].shape
length = raw_data["T"].length

Rd = 287.058 #J/kg/K
Rv = 461.495 #J/kg/K

adjusted_vars = ["U", "V", "W", "Z","P_TOT","RHO_A","RHO_M","T","QVAPOR","REF"]
terrain_vars = ["DHDX", "DHDY"]
adjusted_vars.each {|var| adjusted_data[var] = NArray.float(shape[0],shape[1],shape[2])}
terrain_vars.each {|var| terrain_data[var] = NArray.float(shape[0],shape[1])}

# Writing Time, Lat, Lon, Z, U, V, W, T, Q_vapor, Rho_m and Rho_a to file
puts "Writing to " + fileout + " ..."
puts "Writing to " + fileterrain + " ..."

# puts "Shape = " + shape[0].to_s + shape[1].to_s + shape[2].to_s
bg = File.open(fileout,'w')
terrain = File.open(fileterrain,'w')

(0..shape[0]-1).each do |i|
  (0..shape[1]-1).each do |j|
    lat = raw_data["XLAT"][i, j, timestep]
    lon = raw_data["XLONG"][i, j, timestep]
    terrainHeight = raw_data["HGT"][i, j, timestep]
    # Loop over all dimensions to compute the x-derivative
    # Forward difference at innermost boundary
    if i==0
      terrain_data["DHDX"][i,j] = (-3.0*raw_data["HGT"][i,j,timestep]+4.0*raw_data["HGT"][i+1,j,timestep]-raw_data["HGT"][i+2,j,timestep])/(dx*2)
    # Reverse difference at outermost boundary
    elsif i==shape[0]-1
      terrain_data["DHDX"][i,j] = (3.0*raw_data["HGT"][i,j,timestep] - 4.0*raw_data["HGT"][i-1,j,timestep] + raw_data["HGT"][i-2,j,timestep])/(dx*2)
    # Centered difference at all other grid points
    else
      terrain_data["DHDX"][i,j] = (raw_data["HGT"][i+1,j,timestep] - raw_data["HGT"][i-1,j,timestep])/(dx*2)
    end

    # Forward difference at innermost boundary
    if j==0
      terrain_data["DHDY"][i,j] = (-3.0*raw_data["HGT"][i,j,timestep] + 4.0*raw_data["HGT"][i,j+1,timestep]-raw_data["HGT"][i,j+2,timestep])/(dy*2)
    # Reverse difference at outermost boundary
  elsif j==shape[1]-1
      terrain_data["DHDY"][i,j] = (3.0*raw_data["HGT"][i,j,timestep] - 4.0*raw_data["HGT"][i,j-1,timestep] + raw_data["HGT"][i,j-2,timestep])/(dy*2)
    # Centered difference at all other grid points
    else
      terrain_data["DHDY"][i,j] = (raw_data["HGT"][i,j+1,timestep] - raw_data["HGT"][i,j-1,timestep])/(dy*2)
    end

    terrain.write([lat,lon,terrainHeight,terrain_data["DHDX"][i,j],terrain_data["DHDY"][i,j]].join(" \t")+"\n")

    (0..shape[2]-1).each do |k|


      p_tot = raw_data["P"][i, j, k].to_r+raw_data["PB"][i, j, k].to_r
      e = raw_data["QVAPOR"][i, j, k].to_r*p_tot/0.622
      theta = raw_data["T"][i, j, k].to_r + 300.0
      t = theta*(p_tot/100000.0)**(2/7.0)

      adjusted_data["T"][i, j, k] = t
      adjusted_data["P_TOT"][i, j, k] = p_tot
      adjusted_data["RHO_A"][i, j, k] = (p_tot-e)/(Rd*t)
      adjusted_data["RHO_M"][i, j, k] = (p_tot-e)/(Rd*t)+e/(Rv+t)
      adjusted_data["QVAPOR"][i, j, k] = raw_data["QVAPOR"][i, j, k].to_r*1000.0
      adjusted_data["REF"][i, j, k] = raw_data["REFL_10CM"][i, j, k]

      # Unstagger geopotential heights in vertical
      ph = (raw_data["PH"][i, j, k].to_r + raw_data["PH"][i, j, k+1].to_r)/2.0
      phb = (raw_data["PHB"][i, j, k].to_r + raw_data["PHB"][i, j, k+1].to_r)/2.0
      adjusted_data["Z"][i, j, k] = (ph + phb)/9.81

      # Unstagger U in east-west direction
      adjusted_data["U"][i, j, k] = (raw_data["U"][i, j, k].to_r + raw_data["U"][i+1, j, k].to_r)/2.0

      # Unstagger V in north-south direction
      adjusted_data["V"][i, j, k] = (raw_data["V"][i, j, k].to_r + raw_data["V"][i, j+1, k].to_r)/2.0

      # Unstagger U in east-west direction
      adjusted_data["W"][i, j, k] = (raw_data["W"][i, j, k].to_r + raw_data["W"][i, j, k+1].to_r)/2.0

      # if calc_dist(center_lat,center_lon,lat,lon) <= maxradius then
      if (k == 0) then
        sfcu = raw_data["U10"][i, j, timestep]
        sfcv = raw_data["V10"][i, j, timestep]
        sfcw = 0.0
        sfchgt = 10.0 + terrainHeight
        sfctemp = raw_data["T2"][i, j, timestep]
        sfcpress = raw_data["PSFC"][i, j, timestep]
        sfcq = raw_data["Q2"][i, j, timestep]*1000.0
        sfce = sfcq*sfcpress/622.0
        sfcrhoa = (sfcpress-sfce)/(Rd*sfctemp)
        sfcrhom = sfcrhoa+sfce/(Rv*sfctemp)
        sfcz = adjusted_data["REF"][i, j, 0];
        bg.write([time_input,lat,lon,sfchgt,sfcu,sfcv,sfcw,sfctemp,sfcq,sfcrhoa,sfcrhom,sfcz].join(" \t")+"\n")
      else
        bg.write([time_input,lat,lon,adjusted_data["Z"][i, j, k],\
        adjusted_data["U"][i, j, k],adjusted_data["V"][i, j, k],adjusted_data["W"][i, j, k],\
        adjusted_data["T"][i, j, k],adjusted_data["QVAPOR"][i, j, k],adjusted_data["RHO_A"][i, j, k], adjusted_data["RHO_M"][i, j, k],\
        adjusted_data["REF"][i, j, k]].join(" \t")+"\n")
      end
      # end
    end
  end
end

# # # Clean up
puts "Done :)"
