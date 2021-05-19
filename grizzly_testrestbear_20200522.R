trk.mem = sim.track(mu = 1, kappa = 1, beta.mem = 20, mean.mem = 50, sd.mem = 25, 
                    threshold = 1, n.steps = 500, burnin = 250, gridsize = NULL)
trk.mem = add.resting.point(trk.mem, start = round(runif(1, 250, nrow(trk.mem))), 25)
trk.tofit = sim.controls(trk.mem, n.controls = 100, burnin = 250)
fit = fit.SSF(track = trk.tofit, models = 'memory', export.functions = TRUE,
              init.cond = 3, lb.sdtau = 10)
