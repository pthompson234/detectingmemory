source('sim_functions_20200326.R')
rm(sim.track, sim.and.fit.tracks, make.SSF.raw.data) # clear some bigger functions from memory

require(circular)
require(doParallel)

BEARIDs = list.files('inputdata') # This will work if you have ran grizzly_preapredata_20200423 already
registerDoParallel(cores = length(BEARIDs))

# Note the parallel for loop here. This works just like a for loop if you run on your personal computer and works much more quickly if you run remotely (which I would highly recommend).

results.all = foreach (ID = BEARIDs, .combine = rbind, .errorhandling = 'stop') %dopar% {
  print(paste0('beginning bear ', ID))
  data.raw = read.SSF.raw.data(paste0('inputdata/', ID))
  
  fit = fit.SSF(data = data.raw, init.cond = 3, lb.sdtau = 18, BIC = TRUE,
                mu.stationary = 30, control = list(eval.max = 1000, iter.max = 1000),
                models = c('null', 'resource', 'memory', 'combination'))
  
  cbind(ID = ID, fit)
}

write.csv(results.all, 'results_all_bears.csv')