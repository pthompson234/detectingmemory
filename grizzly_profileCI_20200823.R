source('sim_functions_20200326.R')

require(doParallel)

BEARIDs = list.files('inputdata')
registerDoParallel(cores = length(BEARIDs))

initial.steps.null = rep(1, 4)
initial.steps.res = c(1, 1, 1, 1e-5, 1, 1, 1e-5, 1e-5, 1, 1)
initial.steps.mem = rep(1, 8)
initial.steps.comb = c(1, 1, 1, 1, 1e-5, 1, 1, 1e-5, 1e-5, 1, 1, 1, 1, 1, 1)

dir.create('CI_Info')
sink(file = 'CI_Info/profile_ci.log')

CIs = foreach(ID = BEARIDs, .combine = rbind, .errorhandling = 'remove') %dopar% {
  data.raw = read.SSF.raw.data(paste0('inputs_n50_j30/', ID))
  fit = fit.SSF(data = data.raw, init.cond = 3, lb.sdtau = 18, BIC = TRUE,
                mu.stationary = 30, control = list(eval.max = 1000, iter.max = 1000),
		export.functions = TRUE)
  
  # write.csv(fit, paste0('CI_Info/fit_', ID, '.csv'))
  # If you have the fits already, no need to write them.
  min.BIC.index = which(fit$IC == min(fit$IC))[1]
  
  if (min.BIC.index == 23) {
    #get CI lower and upper bound for each
    CI.comb = t(sapply(X = 1:15, FUN = function(i) {
      prof.lik.CI(fun = L.comb1pk$fn, index = i, tol = initial.steps.comb[i] * 1e-3,
                  optimum = fit$estimate[23:37], gr = L.comb1pk$gr,
                  constraints = matrix(ncol = 2, data = c(fit$const.l[23:37],
                                                          fit$const.u[23:37])),
                  init.step = initial.steps.comb[i], custom.message = paste0(ID, ' comb '), 
                  write = paste0('CI_Info/profiles_20200825/', ID, 'comb', i, '.csv'))
    }))
  } else if (min.BIC.index == 15) {
    #get CI lower and upper bound for each
    CI.mem = t(sapply(X = 1:8, FUN = function(i) {
      prof.lik.CI(fun = L.tau1pk$fn, index = i, tol = initial.steps.mem[i] * 1e-3,
                  optimum = fit$estimate[15:22], gr = L.tau1pk$gr,
                  constraints = matrix(ncol = 2, data = c(fit$const.l[15:22],
                                                          fit$const.u[15:22])),
                  custom.message = paste0(ID, ' mem '), init.step = initial.steps.mem[i],
                  write = paste0('CI_Info/profiles_20200825/', ID, 'mem', i, '.csv'))
      
    }))
  } else if (min.BIC.index == 5) {
    
    #get CI lower and upper bound for each
    CI.res = t(sapply(X = 1:10, FUN = function(i) {
      prof.lik.CI(fun = L.res$fn, index = i, tol = initial.steps.res[i] * 1e-3,
                  optimum = fit$estimate[5:14], gr = L.res$gr,
                  constraints = matrix(ncol = 2, data = c(fit$const.l[5:14],
                                                          fit$const.u[5:14])),
                  init.step = initial.steps.res[i], custom.message = paste0(ID, ' res '), 
                  write = paste0('CI_Info/profiles_20200825/', ID, 'res', i, '.csv'))
    }))
  } else {
    #get CI lower and upper bound for each
    CI.null = t(sapply(X = 1:4, FUN = function(i) {
      prof.lik.CI(fun = L.null$fn, index = i, tol = initial.steps.null[i] * 1e-3,
                  optimum = fit$estimate[1:4], gr = L.null$gr,
                  constraints = matrix(ncol = 2, data = c(fit$const.l[1:4],
                                                          fit$const.u[1:4])),
                  custom.message = paste0(ID, ' null '), init.step = initial.steps.null[i],
                  write = paste0('CI_Info/profiles_20200825/', ID, 'null', i, '.csv'))
    }))
  }
}