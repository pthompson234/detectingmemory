source('sim_functions_20200326.R')

# some things to note here:

# n.sim = 1; this is for ease of simulation. If you would like to run more simulations you are more than welcome to. You will also need to use write.NLM to make more simulated resource landscapes becuase I only provide three here.
# parallel = 0 in sim.and.fit.tracks; if you are choosing to run multiple simulations, you may want to run them remotely and change this setting.
# Use GIF = TRUE for each call to sim.and.fit.tracks if you would like to replicate the generated plots. Unfortunately, you must run one at a time and then move the files so they are not overwritten (this is something I want to fix in the future)

set.seed(123)

n.sim = 1
n.controls = 10
n.steps = 1200 # takes off first 600 steps for burnin and this gets 600 steps of fitting

A = matrix(nrow = 2, data = c(0.85, 0.1, 0.15, 0.9), byrow = TRUE) # Markov state switching matrix

#null model
sim.and.fit.tracks(N = n.sim, gridsize = 100, burnin = 600, n.controls = n.controls, 
                   filename = paste0('results_null_N', n.steps, '_n', n.controls, '.csv'), 
                   n.steps = n.steps, mu = 0.75, kappa = 0.75, beta.mem = 0, beta.res = c(0, 0, 0), 
                   init.cond = 1, res.layers = 'NLMs', BIC = TRUE, markov.matrix = A, parallel = 0, 
                   mu.stationary = 0.01, alpha.exp = 1/60, GIF = TRUE)

# resource-only model
sim.and.fit.tracks(N = n.sim, gridsize = 100, burnin = 600, n.controls = n.controls, 
                   filename = paste0('results_res_N', n.steps, '_n', n.controls, '.csv'), 
                   n.steps = n.steps, mu = 0.75, kappa = 0.75, beta.mem = 0, beta.res = c(7.5, -7.5, 0), 
                   init.cond = 1, res.layers = 'NLMs', BIC = TRUE, markov.matrix = A, parallel = 0, 
                   mu.stationary = 0.01, alpha.exp = 1/60, GIF = TRUE)

#memory-only model
sim.and.fit.tracks(N = n.sim, gridsize = 100, burnin = 600, n.controls = n.controls, 
                   filename = paste0('results_mem_N', n.steps, '_n', n.controls, '.csv'), 
                   n.steps = n.steps, mu = 0.75, kappa = 0.75, beta.mem = 50, mean.mem = 550, 
                   sd.mem = 25,  threshold = 1, beta.res = c(0, 0, 0), init.cond = 1, 
                   res.layers = 'NLMs', BIC = TRUE, markov.matrix = A, parallel = 0, 
                   mu.stationary = 0.01, alpha.exp = 1/60, GIF = TRUE)

#resource-memory model
sim.and.fit.tracks(N = n.sim, gridsize = 100, burnin = 600, n.controls = n.controls, 
                   filename = paste0('results_rm_N', n.steps, '_n', n.controls, '.csv'), 
                   n.steps = n.steps, mu = 0.75, kappa = 0.75, beta.mem = 50, mean.mem = 550, 
                   sd.mem = 25,  threshold = 0, beta.res = c(7.5, -7.5, 0), init.cond = 1, 
                   res.layers = 'NLMs', BIC = TRUE, markov.matrix = A, parallel = 0, 
                   mu.stationary = 0.01, alpha.exp = 1/60, GIF = TRUE)
