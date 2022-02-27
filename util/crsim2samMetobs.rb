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
dr = 3*1000 # in m
dy = 3*1000 # in m

init_year = 2021; init_mon  = 7; init_day  = 30; init_hr   = 12
year = 2021; mon  = 7; day  = 30; hr   = 12

filein = "/Users/tingyu/Dropbox/CSU/Research/PhD/samurai/wrf_20210730/crsim/RF_chill.nc"
path = "/Users/tingyu/Dropbox/CSU/Research/PhD/samurai/wrf_20210730/crsim/"
vars1 = ["DV"]
vars2 = ["Zhh"]
radar = "chill"
radar_lat = 40.44625
radar_lon = -104.63708
radar_Alt = 1.5; # in km
fileout = path+radar+"_crsim.rf"
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
filein1 = path+"RF_"+vars1[0]+"_"+radar+".nc"
filein2 = path+"RF_"+vars2[0]+"_"+radar+".nc"
file1 = NetCDF.open(filein1)
file2 = NetCDF.open(filein2)
puts filein1
puts filein2
raw_data = Hash.new()
dims = ["range","elev","naz"]

dims.each do |dim|
  raw_data[dim] = file1.var(dim)
end

vars1.each do |var|
  raw_data[var] = file1.var(var).get[true,true,true]
end

size = raw_data["DV"].size
rank = raw_data["DV"].rank
shape = raw_data["DV"].shape
length = raw_data["DV"].length
# # puts raw_data["DV"]
# # adjusted_vars = ["VD"]
# # adjusted_vars.each {|var| adjusted_data[var] = NArray.float(shape[0],shape[1],shape[2])}

# Writing range, elev, azimuth, variable to file
puts "Writing to " + fileout + " ..."
nr = file1.var("range").get
nel = file1.var("elev").get
naz = file1.var("azim").get
vd = file1.var(vars1[0]).get[true,true,true]
zhh = file2.var(vars2[0]).get[true,true,true]
# rr = nr.get
puts nr[1].to_s

# # puts "Shape = " + shape[0].to_s + shape[1].to_s + shape[2].to_s
bg = File.open(fileout,'w')
bg.write([time_input, radar_lat, radar_lon, radar_Alt].join(" \t")+"\n")
(0..shape[0]-1).each do |i|
  rr = nr[i].to_s
  (0..shape[1]-1).each do |j|
    el = nel[j].to_s
    (0..shape[2]-1).each do |k|
      az = naz[k].to_s
      vdd = vd[i, j, k].to_s
      zhhh = zhh[i, j, k].to_s
      bg.write([az, el, rr, zhhh, vdd].join(" \t")+"\n")
    end
  end
end

# # # Clean up
puts "Done :)"
