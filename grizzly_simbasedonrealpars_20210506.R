source('sim_functions_20200326.R')

require(raster)
require(circular)
require(TMB)
require(fdrtool)

IDs = c('GF1004', 'GF1008', 'GF1016', 'GF1041',
        'GF1086', 'GF1107', 'GF1130', 'GM1046')
# the 8 bears we did the analysis for

all.pars = read.csv('grizzly_pars_forsim.csv')

for (ID in IDs) {
  pars = as.numeric(all.pars[ID][,1])
  A = matrix(byrow = TRUE, nrow = 2, data = c(pars[13], 1-pars[14], 
                                              1-pars[13], pars[14]))
  # Get the actual initial condition for the bear
  ID.x0 = c(read.csv(paste0('inputdata/', ID, '/x_cases_all.csv'))[1,2],
              read.csv(paste0('inputdata/', ID, '/y_cases_all.csv'))[1,2])
  
  sim.and.fit.tracks(N = 1, gridsize = NULL, burnin = 600, n.res = 6,
                     n.controls = 50, n.steps = 1200, raster = TRUE,
                     filename = paste0('sims_', ID, '.csv'), mu = pars[1],
                     kappa = pars[2], beta.res = pars[3:8], threshold = pars[9],
                     beta.mem = pars[10], mean.mem = pars[11], sd.mem = pars[12],
                     alpha.exp = pars[15], init.cond = 1, same.layers = TRUE,
                     res.layers = paste0('rasters/', ID), BIC = TRUE, x0 = ID.x0,
                     markov.matrix = A, parallel = 25, mu.stationary = 30)
}