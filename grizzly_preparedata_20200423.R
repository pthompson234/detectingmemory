source('sim_functions_20200326.R')

require(raster)

set.seed(123) # We did not set a seed in the analysis, but feel free to do so here for consistency. The data generated for the analysis was generated using this script, and with a different random seed. That being said, fitting the model with this different seed should produce remarkably similar results, a testament to the lack of sensitivity to differing random points.

#Load in rasters
berry = raster('inputs/berries.tif')
fish = raster('inputs/fish.tif')
mammal = raster('inputs/squirrel.tif')
hedysarum = raster('inputs/hedysarum.tif')
town = raster('inputs/towns_m.tif')
cabin = raster('inputs/cabins_m.tif')

resources = stack(berry, fish, mammal, hedysarum, town, cabin)
# If using seasonal resources, we can use this.
# res.bounds = matrix(nrow = 2, data = c(213, 332, 130, 289, 254, 332, 101, 166, 0, 366, 0, 366))

# Here are the 8 bears used in this study (and included in the zip file)
BEARIDs = c('GF1004', 'GF1008', 'GF1016', 'GF1041', 'GF1086', 'GF1107', 'GF1130', 'GM1046')

include = numeric(length(BEARIDs))

#Remove bears that don't have more than a year of data cause we can't use them
for (i in 1:length(BEARIDs)) {
  BEARID = BEARIDs[i]
  
  Z_raw = read.csv(paste0('bears_csv/', BEARID, '.csv'))
  Z_year = as.numeric(substr(Z_raw$Date, 1, 4))
  years = unique(Z_year)
  include[i] = length(years) > 1
}

BEARIDs = BEARIDs[as.logical(include)]

results.all = data.frame()

for (ID in BEARIDs) {
  print(paste0('beginning bear ', ID))
  data.raw = read.csv(paste0('bears_csv/', ID, '.csv'))
  trk.raw = fix.raw.data(data.raw, jitter = 30) # This is necessary to prevent the model from strange outputs (see sim_functions documentation under fix.raw.data for more information)
  #turns it into a track for model fitting
  print('simulating available points')
  n.controls = 50
  trk.tofit = sim.controls(track = trk.raw, n.controls = n.controls, gridsize = NULL, burnin = 2190)
  if (nrow(trk.tofit) - nrow(trk.raw) > n.controls * 50) {
    print('writing raw data to CSV')
    dir.create(paste0('inputdata/', ID))
    make.SSF.raw.data(track = trk.tofit, res = resources, res.bounds = res.bounds,
                      write = paste0('inputdata/', ID), 
                      proj4string = proj4string(resources))
  }
}