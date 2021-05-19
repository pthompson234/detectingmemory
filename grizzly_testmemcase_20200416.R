source('sim_functions_20200326.R')

mu = 11
kappa = 1.1
beta = mu * 30
gridsize = NULL

trk = sim.track(burnin = 875, n.steps = 2500, beta.mem = 0, mean.mem = 750, sd.mem = 30, 
                threshold = 1, mu = mu, kappa = kappa, gridsize = gridsize, n.res = 1,
                print = 50)
trk.tofit = sim.controls(trk, n.controls = 10, gridsize = gridsize, burnin = 875)
fit = fit.SSF(track = trk.tofit, res = NULL, true.values = attr(trk, 'true'), 
              init.cond = 3, models = 'memory', export.functions = TRUE, lb.sdtau = 18)
print(fit)

#we want to see the plot

fun = function(x,y) {L.tau1pk$fn(c(mu, kappa, beta, x, y))}
M = plot.lik.surface(fun, gridsize = 75, lower = c(500, 0.5), upper = c(875, 40),
                     xlab = expression(paste(mu, ' - mean timing of periodicity')),
                     ylab = expression(paste(sigma, ' - variation in periodicity')))


