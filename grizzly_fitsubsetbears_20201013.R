source('sim_functions_20200326.R')
rm(sim.track, sim.and.fit.tracks, make.SSF.raw.data) # clear some bigger functions from memory

require(circular)
require(doParallel)
registerDoParallel(cores = 20)

BEARIDs = list.files('inputdata_subset')

# Note the parallel for loop here. This works just like a for loop if you run on your personal computer and works much more quickly if you run remotely (which I would highly recommend).

fits = foreach(ID = BEARIDs, .combine = rbind, .errorhandling = 'stop') %dopar% {
  data.raw = read.SSF.raw.data(paste0('inputdata_subset/', ID))
  fit = fit.SSF(data = data.raw, init.cond = 3, lb.sdtau = 18, BIC = TRUE,
          mu.stationary = 30, control = list(eval.max = 1000, iter.max = 1000))
  cbind(ID = ID, fit)
}

write.csv(fits, paste0('results_all_subsetted_bears.csv'))