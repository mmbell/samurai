=begin
Created on Thr Apr 18 17:47:18 2013
@author: Annette Foerster (foerster@hawaii.edu)
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
center_lat = 22
center_lon = -80
maxradius = 1000  # in km

init_year = 2005; init_mon  = 9; init_day  = 18; init_hr   = 00
year = 2005; mon  = 9; day  = 20; hr   = 19

filein = "/bora/rita2005/WRF/run2/wrf_d4"
fileout = "./samurai_Background.in"

# ------------------  MODIFY END ------------------

# Convert time to seconds and record as time
init_time = DateTime.new(init_year,init_mon,init_day,init_hr)
time = DateTime.new(year,mon,day,hr)
unixtime = time.to_time.to_i.to_s
timestep = ((time-init_time)*24).to_i 

# Read nc-file
puts "Reading NetCDF for " + time.to_s + " at timestep " + timestep.to_s + " ..."
file = NetCDF.open(filein)

raw_data = Hash.new()
adjusted_data = Hash.new()
dims = ["XTIME","XLAT","XLONG"]
vars = ["U","V","W","T","PH","PHB","P","PB","QVAPOR","REFL_10CM"]
sfcvars = ["U10","V10","T2","PSFC","Q2"]

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

adjusted_vars = ["U", "V", "W", "Z","P_TOT","RHO_A","T","QVAPOR","Z_LINEAR"]
adjusted_vars.each {|var| adjusted_data[var] = NArray.float(shape[0],shape[1],shape[2])}

# Writing Time, Lat, Lon, Z, U, V, W, T, Q_vapor and Rho_a to file
puts "Writing to " + fileout + " ..."

bg = File.open(fileout,'w')

(0..shape[0]-1).each do |i|
  (0..shape[1]-1).each do |j|
    (0..shape[2]-1).each do |k|
      lat = raw_data["XLAT"][i, j, timestep]
      lon = raw_data["XLONG"][i, j, timestep]

      p_tot = raw_data["P"][i, j, k].to_r+raw_data["PB"][i, j, k].to_r
      e = raw_data["QVAPOR"][i, j, k].to_r*p_tot/0.622
      theta = raw_data["T"][i, j, k].to_r + 300.0
      t = theta*(p_tot/100000.0)**(2/7.0)

      adjusted_data["T"][i, j, k] = t
      adjusted_data["P_TOT"][i, j, k] = p_tot
      adjusted_data["RHO_A"][i, j, k] = (p_tot-e)/(287.0*t)
      adjusted_data["QVAPOR"][i, j, k] = raw_data["QVAPOR"][i, j, k].to_r*1000.0
      adjusted_data["Z_LINEAR"][i, j, k] = 10.0**(raw_data["REFL_10CM"][i, j, k]*0.1)

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

      if calc_dist(center_lat,center_lon,lat,lon) <= maxradius then
        if (k == 0) then
          sfcu = raw_data["U10"][i, j, timestep]
          sfcv = raw_data["V10"][i, j, timestep]
          sfcw = 0.0
          sfchgt = 10.0
          sfctemp = raw_data["T2"][i, j, timestep]
          sfcpress = raw_data["PSFC"][i, j, timestep]
          sfcq = raw_data["Q2"][i, j, timestep]*1000.0
          sfce = sfcq*sfcpress/622.0 
          sfcrhoa = (sfcpress-sfce)/(287.0*sfctemp)
          sfcz = adjusted_data["Z_LINEAR"][i, j, 0]; 
          bg.write([unixtime,lat,lon,sfchgt,sfcu,sfcv,sfcw,sfctemp,sfcq,sfcrhoa,sfcz].join("\t")+"\n")
        else
          bg.write([unixtime,lat,lon,adjusted_data["Z"][i, j, k],\
          adjusted_data["U"][i, j, k],adjusted_data["V"][i, j, k],adjusted_data["W"][i, j, k],\
          adjusted_data["T"][i, j, k],adjusted_data["QVAPOR"][i, j, k],adjusted_data["RHO_A"][i, j, k], \
          adjusted_data["Z_LINEAR"][i, j, k]].join("\t")+"\n")
        end
      end
    end
  end
end

# Clean up
file.close
bg.close

puts "Done :)"
