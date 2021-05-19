sim.track = function(mu = 1, kappa = 1, mu.s = 0.01, beta.res = rep(0, n.res), 
                     beta.mem = rep(0, n.mem), mean.mem = rep(1, n.mem), sd.mem = rep(1, n.mem), 
                     alpha.mem = 1, threshold = 1, alpha.exp = 1, markov.matrix = diag(2), 
                     pi.mc = c(0, 1),  gridsize = 100, n.steps = 500, n.random = 1e4, n.res = 3, 
                     n.mem = 1, ls.sparsity = 1, burnin = NULL, x0 = NULL, mem.accuracy = 3, 
                     ns.dist = 'exponential', raster = FALSE, print = 10, res.layers = NULL, 
                     GIF = FALSE, store.pars = TRUE, ...) {
  #Simulate an animal movement track over a random landscape
  #Value - a data.frame (X) with x and y coordinates (as well as path metrics) for each step 
  
  #mu - mean step length for exponential distribution
  #kappa - angular concentration parameter
  #mu.s - mean step length for stationary movements
  #beta.res - selection coefficients for resources
  #beta.mem - selection coefficients for memory peaks
  #mean.mem - average timing of each memory mechanism
  #sd.mem - variation in each memory mechanism
  #alpha.mem - power relationship between memory and selection (NO USE IN THE MANUSCRIPT - IS ALWAYS 1)
  #threshold - default expectation of animal's resources
  #alpha.exp - multiplier for distance transformation
  #markov.matrix - Markov transition matrix for movement states (first component is for stationary)
  #pi.mc - Markov initial condition to sample from (first component is for stationary)
  #gridsize - number of cells in one grid row. If NULL, we don't use a grid.
  #n.steps - number of steps to simulate
  #n.random - number of random steps to choose from at each point. If 0 we use rejection sampling
  #n.res - number of resource layers to make
  #n.mem - number of memory mechanisms (for purposes of this analysis, should always be 1)
  #ls.sparsity - average sparsity of the landscape layers; higher values => worse landscapes
  #burnin - when to start including memory mechanism in simulations
  #x0 - the animal's first point on the grid
  #mem.accuracy - numnber of mem.sd's away from mem.mean to include in Gaussian weighted average
  #raster - FALSE if we are working with matrices, TRUE if we are working with rasters
  #print - integer; prints out the status of the simulation every "print" steps
  #res.layers - NULL if we are to create new layers or a vector of .csv file names
  #GIF - logical; do we save the plot at each time step to make an animated GIF
  #store.pars - logical; should we store the true values for parameters in the attributes
  #... - additional arguments to plot.track
  
  if (mu < 0) stop ('mu must be greater than or equal to 0')
  if (any(kappa < 0)) stop('all components of kappa must be non-negative')
  if (any(mean.mem <= 0)) {
    warning('setting any component of mean.mem to be 0 or less is not recommended')
  }
  if (any(sd.mem <= 0)) {
    warning('setting any component of sd.mem to be 0 or less is not recommended')
  }
  if (mem.accuracy <= 0) stop('mem.accuracy must be greater than 0')
  if (mem.accuracy < 2) {
    warning('setting mem.accuracy to be lower than 2 will heavily impact accuracy')
  }
  if (!is.null(burnin) && burnin > n.steps) {
    warning('setting burnin to be greater than n.steps limits beahvior of simulated animals')
  }
  if (!is.integer(n.random) & n.random < 0) stop('n.random should be a non-negative integer')
  if (length(beta.res) != n.res) stop('beta.res must have length n.res')
  if (length(beta.mem) != n.mem) stop('beta.mem must have length n.mem')
  if (length(mean.mem) != n.mem) stop('mean.mem must have length n.mem')
  if (length(sd.mem) != n.mem) stop('sd.mem must have length n.mem')
  if (norm(dim(markov.matrix) - c(2,2), type = '2') != 0) {
    stop('markov.matrix must be a 2x2 matrix')
  }
  if (norm(colSums(markov.matrix) - c(1,1), type = '2') != 0) {
    stop('markov matrix columns must sum to 1')
  }
  if (raster & !is.null(res.layers)) warning('res.layers must be a list of .tif files or else!')
  if (!raster & !is.null(res.layers)) warning('res.layers must be a list of .csv values or else!')
  
  if (is.null(res.layers) & !is.null(gridsize)) require(NLMR) # random fields
  if (raster) require(raster)
  
  #logical variables detailing what model is currently being simulated
  mod.null = (sum(abs(beta.res)) == 0) & (sum(abs(beta.mem)) == 0)
  mod.res = (sum(abs(beta.res)) != 0) & (sum(abs(beta.mem)) == 0)
  mod.mem = (sum(abs(beta.res)) == 0) & (sum(abs(beta.mem)) != 0)
  #if all of these are false then we are fitting the resource-memory model
  
  if (mod.res | mod.null) threshold = 0 # for true values (otherwise doesn't matter)
  if (!mod.null & !mod.mem & is.null(gridsize) & is.null(res.layers) & !raster) {
    # Important point: if we are simulating from the null or memory-only model, we do not need a gridsize, and the animal can go wherever it wants because there is no computational burden of simulating resources. Otherwise, gridsize needs to be defined.
    stop('Simulations with resource selection require a grid')
  }
  if (GIF & is.null(gridsize) & is.null(res.layers)) {
    warning("You've selected 'GIF = TRUE' but also 'gridsize = NULL' so we can't make a plot")
  }
  if (!is.null(gridsize) | raster) {
    #Make the resource covariate array (B)
    if (raster) {
      #B is a RasterStack
      rasters = list()
      if (is.null(res.layers)) {
        #make random rasters using NLMR
        for (r in 1:n.res) {
          rasters[[r]] = nlm_gaussianfield(gridsize, gridsize)^ls.sparsity
        }
      } else {
        #load in rasters which should be .tif files
        for (r in 1:n.res) {
          rasters[[r]] = raster(res.layers[r])
        }
      }
      B = stack(rasters)
    } else {
      #B is a three-dimensional array
      B = array(dim = c(gridsize, gridsize, n.res))
      if (is.null(res.layers)) {
        #make layers into matrices
        for (r in 1:n.res) {
          B.data = nlm_gaussianfield(gridsize, gridsize)
          B[,,r] = matrix(nrow = gridsize, data = B.data@data@values ^ ls.sparsity)[,gridsize:1]
        }
      } else {
        #load in layers from file
        for (r in 1:n.res) {
          B[,,r] = as.matrix(read.csv(res.layers[r])[,-1])
        }
      }
    }
    B <<- B # B is stored in the global environment for access outside the function
  }
  
  if (!raster) gridsize = dim(B)[1]
  # if res.layers != NULL then we now do have a matrix so let's re-assign gridsize 
  
  X = matrix(nrow = n.steps + 1, ncol = 6, data = 0,
             dimnames = list(0:n.steps, c('t', 'x', 'y', 'theta', 'phi', 'r')))
  #theta - heading; phi - turning angle; r - step length
  
  #initiate the animal at a random point on the grid
  if (is.null(gridsize) & is.null(x0)) {
    X[1, ] = c(0, 0, 0, 0, 0, 0) # this is OK because it's all relative if we don't have a grid
  } else {
    if (!raster & !is.null(x0) && !on.grid(x = x0[1], y = x0[2], gridsize = gridsize)) {
      stop('x0 must be on the grid')
    }
    if (!is.null(x0)) {
      X[1, ] = c(0, x0, 0, 0, 0) 
      # fill path metrics with zeros for the first point (not factored into any analysis)
    } else {
      X[1, ] = c(0, runif(2) * gridsize, 0, 0, 0) # assign point randomly by default
    }
  }
  
  if (!is.null(gridsize) | raster) {
    #A is the linear combination of resource values
    if (!mod.null & !mod.mem) { # if we need A_res make it, otherwise just for plots
      if (raster) {
        A_res.img = A_res = sum(B * beta.res)
      } else {
        A_res.img = A_res = apply(X = B * array(dim = c(gridsize, gridsize, n.res), 
                                                data = rep(beta.res, each = gridsize ^ 2)),
                                  MARGIN = c(1,2), FUN = sum)
      }
    } else {
      if (raster) {
        A_res.img = B[[1]]
        A_res = 0 * A_res.img
      } else {
        A_res.img = B[,,1] # this is really only important for plotting purposes
        A_res = 0 * A_res.img # W(x) == 0 for memory model and we use this for calculations!
      }
    }
  }
  
  #make the Gaussian weights now so as to reduce computational expense
  if (!mod.null & !mod.res) {
    dnorms = matrix(nrow = n.mem, ncol = n.steps)
    for (k in 1:n.mem) {
      dnorms[k, ] = dnorm(1:n.steps, mean.mem[k], sd.mem[k])
    }
  }
  
  #get the GIF-making process ready
  if (GIF) {
    if (!dir.exists('GIFs')) dir.create('GIFs')
    num.digits = ceiling(log10(n.steps)) + 1 # how many digits to index GIF images
    png(paste0('GIFs/z%0', num.digits, 'd.png'))
  }
  
  #define the burnin
  if ((!mod.null & !mod.res) & is.null(burnin)) {
    burnin = max(mean.mem)
  }
  
  pi.mc = pi.mc / sum(pi.mc)
  state = sample(x = c('s', 'ns'), size = 1, prob = pi.mc) # beginning state
  
  for (t in 1:n.steps) {
    if (state == 'ns') {
      if (mod.null | (mod.mem && t < burnin)) {
        #don't need to do any special sampling technique if either a) we are in the null model or b) we are "training" the memory-only model
        if (raster) {
          next.steps = sim.rand.steps(1, mu, kappa, gridsize, X[t, 'x'], X[t, 'y'], X[t, 'theta'],
                                      dist = ns.dist, raster = A_res)
        } else {
          next.steps = sim.rand.steps(1, mu, kappa, gridsize, X[t, 'x'], X[t, 'y'], X[t, 'theta'],
                                      dist = ns.dist)
        }
        X[t+1, ] = c(t, next.steps)
      } else {
        #here we conduct a Monte Carlo simulation at each step, approximating the animal's redistribution kernel by simulating n.random points (simulated from a correlated random walk) and pick based on how attractive (i.e., W) these points are to the animal.
        if (raster) {
          next.steps = sim.rand.steps(n.random, mu, kappa, gridsize, X[t, 'x'], 
                                      X[t, 'y'], X[t, 'theta'], dist = ns.dist,
                                      raster = A_res)
        } else {
          next.steps = sim.rand.steps(n.random, mu, kappa, gridsize, X[t, 'x'], 
                                      X[t, 'y'], X[t, 'theta'], dist = ns.dist)
        }
        if (!mod.res && t >= burnin) {
          #if we are simulating from the memory or resource-memory model and t >= burnin
          W_mem = rep(0, n.random)
          for (k in 1:n.mem) {
            #we judge whether we include a normal weight based on how close it is to mean.mem
            #by excluding faraway values (which are basically 0) we save some computational time
            #it also can't be greater than t obviously
            bools.norm = abs((1:n.steps) - mean.mem[k]) < (sd.mem[k] * mem.accuracy)
            bools.norm = bools.norm & (1:n.steps <= t)
            tau.values = (1:n.steps)[bools.norm] 
            # the values of tau that we will use in the weighted sum
            if (length(tau.values) == 0) stop('No memory to draw from. Try a higher burnin.')
            past.indices = t - tau.values + 1 # the points in the track that we use for the weighted sum
            dnorms.k = dnorms[k, tau.values] / sum(dnorm(tau.values, mean.mem[k], sd.mem[k]))
            if(any(is.infinite(dnorms.k))) {
              stop('Sum of dnorm values is 0. Try a higher sd.mem')
            }
            # Get the animal's distance from each previously visited point on its track (see weighted.dist for more on the function itself)
            dist.k = mapply(weighted.dist, startx = X[past.indices, 'x'], 
                            starty = X[past.indices, 'y'], MoreArgs = list(x = next.steps[,'x'], 
                                                                           y = next.steps[,'y'], 
                                                                           alpha = alpha.mem))
            # Exponential decay with alpha
            dist.k = exp(-alpha.exp * dist.k)
            # Get the resource values from previously visited locations (for resource-memory model)
            if (!is.null(gridsize) | raster) {
              if (raster) {
                past.sp = SpatialPoints(X[past.indices, c('x','y')])
                res.past = extract(A_res, past.sp)
              } else {
                res.past = extract.matrix(X[past.indices, 'x'], X[past.indices, 'y'], A_res)
              }
            } else {
              res.past = 0
            }
            # ###FOR DEBUGGING PURPOSES ONLY
            # if (any(is.na(res.past))) {
            #   print(extract(A_res, SpatialPoints(data.frame(x = X[1, 'x']+0.001, y = X[1, 'y']+0.001))))
            #   stop("stopped")
            # }
            
            dist.k = dist.k %*% ((res.past + threshold) * dnorms.k)
            #quick linear combination of tau portions of the sum
            W_mem = W_mem + beta.mem[k] * dist.k
          }
        } else {
          W_mem = 0
        }
        if (!mod.mem) {
          # we have to add on resources
          if (raster) {
            next.steps.sp = SpatialPoints(next.steps[, c('x','y')])
            W_res = extract(A_res, next.steps.sp)
          } else {
            W_res = extract.matrix(next.steps[, 'x'], next.steps[, 'y'], A_res)
          }
        } else {
          W_res = 0
        }
        W_all = exp(W_res + W_mem)
        if (all(is.na(W_all))) stop('every value of W_all is NA')
        W_all[is.na(W_all)] = 0 # important for rasters
        # ### FOR DEBUGGING PURPOSES ONLY
        # if (sum(W_all) == 0) {
        #   write.csv(W_res, 'W_res.csv')
        #   write.csv(W_mem, 'W_mem.csv')
        # }
        if (sum(W_all) == 0) stop('empty probability vector!')
        W_all_prob = W_all / sum(W_all)
        step = sample(1:n.random, size = 1, prob = W_all_prob) 
        # sample randomly based on probabilities
        X[t+1, ] = c(t, next.steps[step, ])
      }
    } else {
      #don't need to do any special sampling technique in the stationary state
      if (raster) {
        next.steps = sim.rand.steps(1, mu.s, 0, gridsize, X[t, 'x'], X[t, 'y'], X[t, 'theta'],
                                    dist = 'gaussian', raster = A_res)
      } else {
        next.steps = sim.rand.steps(1, mu.s, 0, gridsize, X[t, 'x'], X[t, 'y'], X[t, 'theta'],
                                    dist = 'gaussian')
      }
      X[t+1, ] = c(t, next.steps)
    }
    
    if (!is.null(gridsize) & GIF) {
      #make the image for this point in time
      plot.track(track = X[1:(t+1), c('x', 'y')] / gridsize, bg = A_res.img, 
                 pt.col = 'red', main = paste0('t = ', t), ...)
    }
    
    #get next state
    current.state.vector = as.numeric(c(state == 's', state == 'ns'))
    next.state.probs = markov.matrix %*% current.state.vector
    state = sample(x = c('s', 'ns'), size = 1, prob = next.state.probs)
    
    if (print > 0 && t %% print == 0) print(paste0('Completed step ', t))
  }
  
  if (GIF) dev.off()
  
  if (store.pars) {
    pars.list = list()
    pars.list[['null']] = c(mu, kappa, log(markov.matrix[1,1] / (1-markov.matrix[1,1])),
                            log(markov.matrix[2,2] / (1-markov.matrix[2,2])))
    pars.list[['resource']] = c(mu, kappa, beta.res, 
                                log(markov.matrix[1,1] / (1-markov.matrix[1,1])),
                                log(markov.matrix[2,2] / (1-markov.matrix[2,2])))
    pars.list[['memory (1-peak)']] = c(mu, kappa, beta.mem[1], mean.mem[1], sd.mem[1],
                                       log(markov.matrix[1,1] / (1-markov.matrix[1,1])),
                                       log(markov.matrix[2,2] / (1-markov.matrix[2,2])),
                                       log(alpha.exp))
    pars.list[['memory (n-peak)']] = c(mu, kappa, beta.mem, mean.mem, sd.mem,
                                       log(markov.matrix[1,1] / (1-markov.matrix[1,1])),
                                       log(markov.matrix[2,2] / (1-markov.matrix[2,2])),
                                       log(alpha.exp))
    pars.list[['memory (power)']] = c(mu, kappa, beta.mem, mean.mem, sd.mem, alpha.mem,
                                      log(markov.matrix[1,1] / (1-markov.matrix[1,1])),
                                      log(markov.matrix[2,2] / (1-markov.matrix[2,2])))
    pars.list[['combination (1-peak)']] = c(mu, kappa, threshold, beta.res, beta.mem[1], 
                                            mean.mem[1], sd.mem[1],
                                            log(markov.matrix[1,1] / (1-markov.matrix[1,1])),
                                            log(markov.matrix[2,2] / (1-markov.matrix[2,2])),
                                            log(alpha.exp))
    pars.list[['combination (n-peak)']] = c(mu, kappa, threshold, beta.res, beta.mem, 
                                            mean.mem, sd.mem,
                                            log(markov.matrix[1,1] / (1-markov.matrix[1,1])),
                                            log(markov.matrix[2,2] / (1-markov.matrix[2,2])),
                                            log(alpha.exp))
    pars.list[['combination (power)']] = c(mu, kappa, threshold, beta.res, beta.mem, 
                                           mean.mem, sd.mem, alpha.mem,
                                           log(markov.matrix[1,1] / (1-markov.matrix[1,1])),
                                           log(markov.matrix[2,2] / (1-markov.matrix[2,2])))
  }
  
  attr(X, 'true') = pars.list
  attr(X, 'gridsize') = gridsize
  X
}

sim.rand.steps = function(n, mu, kappa, gridsize, startx, starty, start.theta,
                          dist = 'exponential', raster = NULL) {
  #n - number of points to simulate
  #mu - mean step length (no matter what distribution this parameter indicates the mean)
  #kappa - angular concentration parameter
  #gridsize - number of grid cells (points cannot fall off the grid); can also be NULL
  #startx, starty - coordinates of original point
  #start.theta - bearing the animal arrived at (startx, starty) on
  #dist - step length distribution to be used - either 'exponential' or 'gaussian' or 'rayleigh'
  #raster - if not NULL, a RasterStack or RasterLayer
  require(circular) #rvonmises
  
  if (!on.grid(startx, starty, gridsize, raster)) stop('animal is not starting on the grid')
  
  zero.circular = circular(0) # by doing this, R will not convert every single number to a circular (this is very time consuming)
  
  if (dist == 'exponential') {
    dist.fn = rexp
    dist.par = 1/mu
  } else if (dist == 'gaussian') {
    require(fdrtool) #rhalfnorm
    dist.fn = rhalfnorm
    dist.par = sqrt(pi/2) / mu
  } else if (dist == 'rayleigh') {
    require(VGAM) #rrayleigh
    dist.fn = rrayleigh
    dist.par = mu / sqrt(pi/2)
  } else {
    stop('invalid distribution parameter')
  }
  
  next.steps = matrix(nrow = n, ncol = 5)
  next.steps[,5] = dist.fn(n, dist.par)
  next.steps[,4] = rvonmises(n, zero.circular, kappa) * sign(next.steps[,5])
  next.steps[,3] = (next.steps[,4] + start.theta) %% (2*pi)
  next.steps[,1] = startx + next.steps[,5] * cos(next.steps[,3])
  next.steps[,2] = starty + next.steps[,5] * sin(next.steps[,3])
  while(any(!on.grid(next.steps[,1], next.steps[,2], gridsize, raster))) {
    # We are not using reflective boundaries, but we prevent the simulated animal from leaving the grid by re-simulating any step until it is on the grid. As a result, this produces behavior that confines the animal within the grid, although it does not "bounce" off as one might expect from reflective boundaries.
    indices = which(!on.grid(next.steps[,1], next.steps[,2], gridsize, raster))
    n.indices = length(indices)
    next.steps[indices,5] = dist.fn(n.indices, dist.par)
    next.steps[indices,4] = rvonmises(n.indices, zero.circular, kappa)
    next.steps[indices,3] = (next.steps[indices,4] + start.theta) %% (2*pi)
    next.steps[indices,1] = startx + next.steps[indices,5] * cos(next.steps[indices,3])
    next.steps[indices,2] = starty + next.steps[indices,5] * sin(next.steps[indices,3])
  }
  
  colnames(next.steps) = c('x', 'y', 'theta', 'phi', 'r')
  
  next.steps
}

on.grid = function(x, y, gridsize, raster = NULL) {
  #x, y - coordinates
  #gridsize - maximum extent of x and y - can also be NULL
  #raster - NULL or a RasterLayer
  
  #returns a logical describing whether the point (x,y) is on a grid of size gridsize
  
  if (!is.null(raster)) {
    return(!is.na(extract(raster, SpatialPoints(data.frame(x = x, y = y)))))
  } else {
    min.x = min.y = 0
    max.x = max.y = gridsize
  }
  
  if (!is.null(gridsize) & !is.null(raster)) stop('gridsize or raster should be NULL')
  if (length(x) != length(y)) stop('x and y must have the same length')
  if (is.null(gridsize) & is.null(raster)) return(rep(TRUE, length(x)))
  return(x > min.x & y > min.y & x < max.x & y < max.y)
}

extract.matrix = function(x, y, matrix) {
  #x, y - coordinates
  #matrix - square matrix where x and y match as coordinates to matrix indices
  
  # Akin to the function raster::extract but for a matrix representing a grid. All values round up, i.e., the point (0.2, 0.7) will give the matrix index [1,1] (of course noting that R uses one-based indexing)
  
  if (nrow(matrix) != ncol(matrix)) stop('matrix must be square')
  if (any(!on.grid(x, y, nrow(matrix)))) stop('some points are off the matrix grid')
  
  x = ceiling(x)
  y = ceiling(y) #have to make both of these integers so they can be indexed
  unlist(lapply(1:length(x), function(i) {
    matrix[x[i], y[i]]
  }))
}

weighted.dist = function(x, y, startx, starty, alpha = 1, beta = 0) {
  #x, y - vector of current points
  #startx, starty - coordinates for start point (1 point)
  #alpha - power of distance (1 implies traditional Euclidean distance)
  #beta - shift in distance (useful for alpha < 0); 0 implies traditional distance
  
  # This is sort of deprecated but was generated for the testing of other model forms. As of now it returns basic Euclidean distance.
  if (length(x) != length(y)) stop('x and y should have the same length')
  
  (beta + sqrt((x - startx)^2 + (y - starty)^2))^alpha
}

plot.track = function(track, bg, bg.col = hcl.colors(64, 'Ag_GrnYl'), 
                      pt.col = 'black', ln.col = 'black', last.col = pt.col, 
                      seasons.col = NULL, ...) {
  # Plots a simulated GPS track of animal movements
  
  #track - data frame with 2 columns (x, y)
  #bg - background matrix to be plotted using image.plot
  #bg.col - vector of colors to be included in image.plot
  #pt.col - the color of the animal movement points (one color)
  #last.col - the color of the last point (default is just the same color, but if one would like the last point to stand out this is here)
  #seasons.col - NULL or a range of colors for the time of period (one for each point)
  #... - additional arguments to image.plot
  require(fields) #image.plot
  
  if ('RasterLayer' %in% class(bg)) {
    #we have to convert to matrix to use image.plot
    data.bg = rev(bg@data@values)
    bg = matrix(nrow = sqrt(length(data.bg)), data = data.bg, byrow = TRUE)
  }
  
  image.plot(bg, col = bg.col, ...)
  if (!is.null(seasons.col)) {
    points(track[-nrow(track), 1], track[-nrow(track), 2], pch = 16, 
           col = seasons.col[-length(seasons.col)])
    points(track[nrow(track), 1], track[nrow(track), 2], pch = 16, 
           col = seasons.col[length(seasons.col)])
  } else {
    points(track[-nrow(track), 1], track[-nrow(track), 2], pch = 16, col = pt.col)
    points(track[nrow(track), 1], track[nrow(track), 2], pch = 16, col = last.col)
  }
  lines(track[, 1], track[, 2], col = ln.col)
  
}

sim.controls = function(track, n.controls = 10, gridsize = attr(track, 'gridsize'),
                        burnin = 1, move.pars = NULL, raster = NULL) {
  #track - matrix of the form from sim.track (6 columns) - whether you intend to fit them or not
  #n.controls - number of simulated points per observed point
  #gridsize - number of cells in one grid row; control steps cannot fall off the grid
  #burnin - value of 't' for which we begin fitting points (i.e., simulating controls)
  #move.pars - either NULL or a parameter vector to simulate movement tracks from
  #raster - do we have a raster? if not NULL, then a RasterStack or RasterLayer
  
  #Return value - a data.frame with the observed (or "case") and simulated ("available" or "control") points
  if (!is.null(gridsize) & !is.null(raster)) stop('gridsize or raster should be NULL')
  
  require(circular)
  track.withcontrol = cbind(track, obs = 1) # 'obs' column is 1 for case, 0 for control
  
  # which points do we actually fit? See more about codes in fix.raw.data
  if ("code" %in% colnames(track)) {
    fit.indices = (track[, 'code'] == 1 | track[, 'code'] == 2) & track[,'t'] >= burnin
  } else {
    fit.indices = track[,'t'] >= burnin 
  }
  if (is.null(move.pars)) {
    mean.steplength = mean(track[fit.indices, 'r'])
    kappa.turningangle = mle.vonmises(track[fit.indices, 'phi'])$kappa
  } else {
    mean.steplength = move.pars[1]
    kappa.turningangle = move.pars[2]
  }
  
  for (i in (1:nrow(track))[fit.indices]) {
    # Simulate random steps based on the correlated random walk distribution
    next.steps = sim.rand.steps(n = n.controls, mu = mean.steplength, kappa = kappa.turningangle,
                                gridsize = gridsize, startx = track[i-1, 'x'], raster = raster,
                                starty = track[i-1, 'y'], start.theta = track[i-1, 'theta'])
    if ("code" %in% colnames(track)) next.steps = cbind(next.steps, 
                                                        track[i, c('code', 'datetime')])
    # Add to the data frame
    track.withcontrol = rbind(track.withcontrol, cbind(t = track[i, 't'], next.steps, obs = 0))
  }
  
  #get indices of cases we didn't fit - all cases that are FALSE for fit.indices and not controls
  cases.nofit.indices = c(!fit.indices, rep(FALSE, nrow(track.withcontrol) - length(fit.indices)))
  track.withcontrol[cases.nofit.indices, 'obs'] = 2
  
  attr(track.withcontrol, 'true') = attr(track, 'true')
  track.withcontrol
}

add.resting.point = function(track, start, length, sd = 0, 
                             gridsize = attr(track, 'gridsize')) {
  # Used to test the effect of long-term stationary periods on simulated data (not used in any of the analysis here, now that the model has a hidden Markov component).
  
  #track - raw movement data (i.e., returned by sim.track) to add resting point to
  #start - at what time index do we start the resting point
  #length - how long does the animal sit there for
  #sd - degree of stochasticity (~GPS error) in resting locations. Should be small.
  #gridsize - size of grid
  
  if (length == 0) return(track)
  
  rest.steps = sim.track(mu = sd, kappa = 0, gridsize = gridsize, n.steps = length,
                         x0 = track[start, c('x', 'y')], print = FALSE)
  #update the first step after resting
  if (start < nrow(track)) {
    new.theta.after.rest = atan2(track[start+1, 'y']-rest.steps[nrow(rest.steps), 'y'],
                                 track[start+1, 'x']-rest.steps[nrow(rest.steps), 'x'])
    new.theta.after.rest = new.theta.after.rest %% (2*pi)
    new.phi.after.rest = new.theta.after.rest - rest.steps[nrow(rest.steps), 'theta']
    new.r.after.rest = norm(track[start+1, c('x','y')] - rest.steps[nrow(rest.steps), c('x','y')],
                            type = '2')
    track[start + 1, c('theta', 'phi', 'r')] = c(new.theta.after.rest, new.phi.after.rest,
                                                 new.r.after.rest)
  }
  
  #update the first step within the resting
  rest.steps[1, 'phi'] = (- track[start, 'phi']) %% (2*pi)
  
  track.withrest = rbind(track[1:start, ], rest.steps, track[(start+1):nrow(track), ])
  track.withrest[,'t'] = 0:(nrow(track.withrest) - 1)
  attr(track.withrest, 'true') = attr(track, 'true')
  attr(track.withrest, 'gridsize') = attr(track, 'gridsize')
  
  track.withrest
}

fit.SSF = function(track = NULL, res = NULL, true.values = attr(track, 'true'), n.peaks = 1, 
                   power = FALSE, init.cond = 1, data = NULL, 
                   models = c('null', 'resource', 'memory', 'combination'),
                   proj4string = CRS(as.character(NA)), export.functions = FALSE, 
                   lb.sdtau = 1, BIC = FALSE, init.include = FALSE, rm.s = FALSE,
                   fixed.par = NULL, fp.comb = NULL, mu.stationary = 0.01, 
                   delta = NULL, b0 = 1, ...) {
  # The big function. This fits the model to animal movement data (track) and produces outputs. Makes calls to C++ via Template Model Builder.
  
  #track - data (including observed and simulated available steps), including cases that are not being fit (6 columns)
  #res - list of matrices (or a RasterStack) for each environmental layer to be fit
  #true.values - the true values for the models (a list of each parameter name as in sim.track)
  #n.peaks - expected no. of peaks - dictates how many memory models are fit (ALWAYS 1 IN THIS ANALYSIS)
  #power - logical dictating whether the memory and combination power models are fit (ALWAYS FALSE IN THIS ANALYSIS)
  #init.cond - for memory models, the number of initial optimization guesses for memory parameters
  #data - NULL unless the user provides model data so it does not need to be extracted (With real data, this is usually how we pass everything in)
  #models - which models do we fit? Should include some combination of 'null', 'resource', 'memory', or 'combination'. Note that 'combination' is equivalent to the resource-memory model here.
  #proj4string - if there are rasters, use this coordinate projection
  #export.functions - logical; should we export likelihood functions to the global environment
  #lb.sdtau - lower bound for sigma in optimization routine
  #BIC - computes BIC for each model if TRUE, AIC otherwise
  #init.include - include model fits for all initial conditions or just the best one? Meaningless if init.cond == 1
  #rm.s - do we remove stationary points from memory calculation? Always FALSE in this analysis (experimental)
  #fixed.par - do we fix a parameter in memory to profile likelihood? Used in calculation of confidence intervals only.
  #fp.comb - fix resource parameter for combination model (always NULL in this analysis - used for testing)
  #mu.stationary - step length for quasi stationary state (parameter is fixed in estimation)
  #delta - starting state probability for each state
  #b0 - threshold parameter value in resource-memory model if fixed (if it is not fixed then it means nothing)
  #... - additional arguments to nlminb calls (e.g., control)
  
  #Return value: data.frame with parameter estimates and AIC / BIC
  require(TMB) # fitting models
  require(circular) # mle.vonmises
  
  if (('resource' %in% models | 'combination' %in% models) & (is.null(res) & is.null(data))) {
    stop("if model includes resource selection then you must include a res!")
  }
  if (!any('null' %in% models, 'resource' %in% models, 
           'memory' %in% models, 'combination' %in% models)) {
    stop('must include null, resource, memory, or combination in "models"')
  }
  
  results = data.frame()
  
  if (is.null(data)) {
    case.indices = which(track[, 'obs'] == 1)
    control.indices = which(track[, 'obs'] == 0)
    
    r.cases = track[case.indices, 'r']
    strata.controls = track[control.indices, 't']
  } else {
    #we just have to load the data in from the list. We need these names otherwise it won't work!
    r.cases = data$r.cases
    strata.controls = data$strata.controls
  }
  
  if (is.null(true.values)) {
    true.values = list('null' = 0, 'resource' = 0, 'memory (1-peak)' = 0,
                       'memory (n-peak)' = 0, 'memory (power)' = 0,
                       'combination (1-peak)' = 0, 'combination (n-peak)' = 0,
                       'combination (power)' = 0)
  }
  
  burnin = min(strata.controls) # hypothetically at what time did we start fitting the data
  fixed.par.value = fixed.par
  fixed.par = !is.null(fixed.par)
  fp.comb.value = fp.comb
  fp.comb = !is.null(fp.comb)
  
  #this code must be ran in a directory where these files exist already!
  if ('null' %in% models) {
    # Generates the likelihood function for the null model based on input data
    L.null = generateFunction(track = track, res = res, true.values = true.values, data = data,
                              model = 'null', proj4string = proj4string, fixed.par = fixed.par.value,
                              mu.stationary = mu.stationary, delta = delta)
    if (export.functions) L.null <<- L.null
    # Fit the null model to the data
    mle.null = nlminb(L.null$par, L.null$fn, L.null$gr, lower = c(0, 0, -10, -10),
                      upper = c(Inf, Inf, Inf, Inf), ...)
    # Calculate information criteria
    IC.null = ifelse(BIC, 2*mle.null$objective + log(length(r.cases))*length(unlist(L.null$par)), 
                     2*mle.null$objective + 2*length(unlist(L.null$par)))
    # Compile results into a data frame
    results.null = data.frame(model = 'null', IC = IC.null, true = true.values[['null']], 
                              estimate = mle.null$par, se = sdreport(L.null)$sd, 
                              msg = mle.null$message, const.l = c(0, 0, -10, -10), const.u = c(Inf, Inf, Inf, Inf))
    if (!export.functions) rm(L.null)
    results = rbind(results, results.null)
  }
  
  if ('resource' %in% models) {
    # These if statements are always quite similar (note that they are not mutually exclusive - one can fit all four models in one function call) - but they get more complex as the models do
    L.res = generateFunction(track = track, res = res, true.values = true.values, data = data,
                             model = 'resource', proj4string = proj4string, fixed.par = fixed.par.value,
                             mu.stationary = mu.stationary, delta = delta)
    n.res = length(L.res$par)-4
    if (export.functions) L.res <<- L.res
    if (n.res == 6) {
      #we're working with the real bear data, which has large values for distances that require more stringent bounds (they are all in meters in the raw data; keep in mind estimates for these parameters never even come close to 0.1)
      lower.res = c(0, 0, -Inf, -10, -Inf, -Inf, -10, -10, -10, -10)
      upper.res = c(Inf, Inf, Inf, 10, Inf, Inf, 10, 10, 10, 10)
    } else {
      lower.res = c(0,0,rep(-Inf, n.res), rep(-10, 2))
      upper.res = c(rep(Inf, n.res+2), rep(10, 2))
    }
    mle.res = nlminb(L.res$par, L.res$fn, L.res$gr, 
                     lower = lower.res, upper = upper.res, ...)
    IC.res = ifelse(BIC, 2*mle.res$objective + log(length(r.cases))*length(unlist(L.res$par)), 
                    2*mle.res$objective + 2*length(unlist(L.res$par)))
    results.res = data.frame(model = 'resource', IC = IC.res, true = true.values[['resource']], 
                             estimate = mle.res$par, se = sdreport(L.res)$sd, 
                             msg = mle.res$message, const.l = lower.res,
                             const.u = upper.res)
    if (!export.functions) rm(L.res)
    results = rbind(results, results.res)
  }
  
  if ('memory' %in% models | 'combination' %in% models) {
    init.cond.x = seq(burnin, 5, length.out = init.cond)
    #I say 5 as a bit of an arbitrary choice here - trying to get rid of cases where it's 1
    b.mem.low = 0 # may change to -Inf
  }
  
  if ('memory' %in% models) {
    L.tau1pk = generateFunction(track = track, res = res, true.values = true.values, data = data,
                                model = 'memory', proj4string = proj4string, fixed.par = fixed.par.value,
                                mu.stationary = mu.stationary, delta = delta, rm.s = rm.s)
    if (export.functions) L.tau1pk <<- L.tau1pk
    if (fixed.par) {
      # ONLY USED FOR TESTING
      lower.tau1pk = c(0,0,5,lb.sdtau, -10, -10, -Inf)
      upper.tau1pk = c(Inf,Inf,burnin,burnin, Inf, Inf, log(1e10))
    } else {
      #include beta in lower bounds
      mean.steplength = L.tau1pk$par[1]
      # See discussion about this in the MS
      upper.alpha = log(1/mean.steplength)
      # Note the difference between mean.steplength, which incorporates stationary steps, and the step length parameter we estimate
      lower.tau1pk = c(0, 0, b.mem.low, 5, lb.sdtau, -10, -10, -Inf)
      upper.tau1pk = c(Inf, Inf, 1000, burnin, Inf, Inf, Inf, upper.alpha)
    }
    # Since we may try multiple initial conditions for the memory models, we must get a list of models and compare them all, picking the one with the best likelihood
    mle.tau1pk = lapply(X = 1:(init.cond), FUN = function(x) {
      nlminb(c(L.tau1pk$par[1:2], rep(50, !fixed.par), init.cond.x[x], lb.sdtau * 2, 0, 0, upper.alpha), 
             L.tau1pk$fn, L.tau1pk$gr, lower = lower.tau1pk, upper = upper.tau1pk, ...)
    })
    logLiks.tau1pk = unlist(lapply(X = mle.tau1pk, FUN = function(m) {m$objective}))
    if (init.include) {
      if (BIC) {
        ICs.tau1pk = 2*logLiks.tau1pk + log(length(r.cases))*length(unlist(L.tau1pk$par))
      } else {
        ICs.tau1pk = 2*logLiks.tau1pk + 2*length(unlist(L.tau1pk$par))
      }
      results.tau1pk = data.frame(model = 'memory (1-peak)', 
                                  IC = rep(ICs.tau1pk, each = length(L.tau1pk$par)), 
                                  true = rep(true.values[['memory (1-peak)']], each = init.cond^2),
                                  estimate = unlist(lapply(X = mle.tau1pk, 
                                                           FUN = function(m) {m$par})),
                                  se = unlist(lapply(X = mle.tau1pk, 
                                                     FUN = function(m) {
                                                       sdreport(L.tau1pk, par.fixed = m$par)$sd
                                                     })),
                                  msg = rep(unlist(lapply(X = mle.tau1pk, 
                                                          FUN = function(m) {m$message})), 
                                            each = length(L.tau1pk$par)))
    } else {
      #find the intial condition that produced the best results and go with those parameters
      if (any(logLiks.tau1pk > length(r.cases))) {
        #do [1] in case there's a tie, don't want a list of 2 here
        # Note we use min - the TMB functions return negative log likelihood!
        mle.tau1pk = mle.tau1pk[[which(logLiks.tau1pk == min(logLiks.tau1pk[logLiks.tau1pk > length(r.cases)]))[1]]]
      } else {
        #if none of the fits are good just pick the "best" bad one
        mle.tau1pk = mle.tau1pk[[which(logLiks.tau1pk == min(logLiks.tau1pk))[1]]]
      }
      IC.tau1pk = ifelse(BIC, 2*mle.tau1pk$objective + log(length(r.cases))*length(unlist(L.tau1pk$par)), 
                         2*mle.tau1pk$objective + 2*length(unlist(L.tau1pk$par)))
      results.tau1pk = data.frame(model = 'memory (1-peak)', IC = IC.tau1pk, true = true.values[['memory (1-peak)']], 
                                  estimate = mle.tau1pk$par, se = sdreport(L.tau1pk, par.fixed = mle.tau1pk$par)$sd, 
                                  msg = mle.tau1pk$message, const.l = lower.tau1pk, const.u = upper.tau1pk)
    }
    if (!export.functions) rm(L.tau1pk)
    results = rbind(results, results.tau1pk)
    
    if (n.peaks > 1) {
      # Do nothing!
    }
    
    if (power) {
      # Do nothing!
    }
  }
  
  if ('combination' %in% models) {
    L.comb1pk = generateFunction(track = track, res = res, true.values = true.values, data = data,
                                 model = 'combination', proj4string = proj4string, 
                                 fixed.par = fixed.par.value, mu.stationary = mu.stationary, 
                                 delta = delta, b0 = b0)
    
    n.res = length(L.comb1pk$par) - 9
    
    if (export.functions) L.comb1pk <<- L.comb1pk
    if (fixed.par) {
      lower.comb1pk = c(0,0,-Inf, rep(-Inf, n.res-1), b.mem.low, 5,lb.sdtau, -10, -10, -Inf)
      upper.comb1pk = c(Inf, Inf, Inf, rep(Inf, n.res-1), Inf, burnin, burnin, Inf, Inf, log(1e10))
    } else {
      mean.steplength = L.comb1pk$par[1]
      upper.alpha = log(1/mean.steplength)
      
      if (n.res == 6) {
        #working with bear data
        lower.comb1pk = c(0, 0, -Inf, -Inf, -10, -Inf, -Inf, -10, -10, 
                          b.mem.low, 5, lb.sdtau, -10, -10, -Inf)
        upper.comb1pk = c(Inf, Inf, Inf, Inf, 10, Inf, Inf, 10, 10, 
                          1000, burnin, burnin, Inf, Inf, upper.alpha)
        
      } else {
        lower.comb1pk = c(0,0,rep(-Inf, n.res+1), b.mem.low, 1, lb.sdtau, -10, -10, -Inf)
        upper.comb1pk = c(rep(Inf, n.res+3), 1000, burnin, burnin, Inf, Inf, upper.alpha)
      }
    }
    mle.comb1pk = lapply(X = 1:(init.cond), FUN = function(x) {
      start.comb1pk = c(L.comb1pk$par[1:(3+n.res-fp.comb)], 50, init.cond.x[x], 
                        lb.sdtau * 2, 0, 0, min(0, upper.alpha))
      
      # FOR DEBUGGING PURPOSES ONLY
      # print(paste0('Length of starting parameter: ', length(start.comb1pk))
      # print(paste0('Length of function: ', length(L.comb1pk$par)))
      
      nlminb(start.comb1pk, L.comb1pk$fn, L.comb1pk$gr, lower = lower.comb1pk, 
             upper = upper.comb1pk, ...)
    })
    logLiks.comb1pk = unlist(lapply(X = mle.comb1pk, FUN = function(m) {m$objective}))
    if (init.include) {
      if (BIC) {
        ICs.comb1pk = 2*logLiks.comb1pk + log(length(r.cases))*length(L.comb1pk$par)
      } else {
        ICs.comb1pk = 2*logLiks.comb1pk + 2*length(L.comb1pk$par)
      }
      results.comb1pk = data.frame(model = 'combination (1-peak)', IC = rep(ICs.comb1pk, 
                                                                            each = length(L.comb1pk$par)), 
                                   true = rep(true.values[['combination (1-peak)']], each = init.cond^2),
                                   estimate = unlist(lapply(X = mle.comb1pk, FUN = function(m) {m$par})),
                                   se = unlist(lapply(X = mle.comb1pk, FUN = function(m) {
                                     sdreport(L.comb1pk, par.fixed = m$par)$sd
                                   })), 
                                   msg = rep(unlist(lapply(X = mle.comb1pk, FUN = function(m) {m$message})), 
                                             each = length(L.comb1pk$par)))
    } else {
      #find the intial condition that produced the best results and go with those parameters
      if (any(logLiks.comb1pk > length(r.cases))) {
        #do [1] in case there's a tie, don't want a list of 2 here
        mle.comb1pk = mle.comb1pk[[which(logLiks.comb1pk == min(logLiks.comb1pk[logLiks.comb1pk > length(r.cases)]))[1]]]
      } else {
        #if none of the fits are good just pick the "best" bad one
        mle.comb1pk = mle.comb1pk[[which(logLiks.comb1pk == min(logLiks.comb1pk))[1]]]
      }
      IC.comb1pk = ifelse(BIC, 2*mle.comb1pk$objective + log(length(r.cases))*length(L.comb1pk$par), 
                          2*mle.comb1pk$objective + 2*length(L.comb1pk$par))
      results.comb1pk = data.frame(model = 'combination (1-peak)', IC = IC.comb1pk, 
                                   true = true.values[['combination (1-peak)']], 
                                   estimate = mle.comb1pk$par, 
                                   se = sdreport(L.comb1pk, par.fixed = mle.comb1pk$par)$sd, 
                                   msg = mle.comb1pk$message, const.l = lower.comb1pk, 
                                   const.u = upper.comb1pk)
    }      
    if (!export.functions) rm(L.comb1pk)
    results = rbind(results, results.comb1pk)
    
    if (n.peaks > 1) {
      # Do nothing!
    }
    
    if (power) {
      # Do nothing!
    }
  }
  
  if (is.null(true.values)) results = subset(results, select = -c(true))
  cbind(results, CL = results$estimate - 2*results$se, CU = results$estimate + 2*results$se)
  
}


generateFunction = function(track = NULL, res = NULL, true.values = attr(track, 'true'), data = NULL, 
                            model = 'null', proj4string = CRS(as.character(NA)), rm.s = FALSE,
                            fixed.par = NULL, fp.comb = NULL, mu.stationary = 0.01, delta = NULL,
                            b0 = 0) {
  # Given a model and some data, we generate the likelihood function and gradient for the appropriate model.
  
  #track - case and controls for data, including cases that are not being fit (6 columns)
  #res - list of matrices (or a RasterStack) for each environmental layer to be fit
  #true.values - the true values for the models (a list of each parameter name as in sim.track)
  #data - NULL unless the user provides model data so it does not need to be extracted
  #model - which models do we want the function for?
  #proj4string - if there are rasters, use this projection
  #rm.s - do we remove stationary points from memory calculation? Always FALSE in this analysis (experimental)
  #fixed.par - do we fix a parameter in memory to profile likelihood?
  #mu.stationary - step length for quasi stationary state (parameter is fixed in estimation)
  #delta - starting state probability for each state - (0,1) is the old model
  #b0 - value for threshold parameter if it's fixed in the model (means nothing in this case since we are not using this)
  
  #Return value: function returned by MakeADFun
  
  require(TMB) # fitting models
  require(circular) # mle.vonmises
  
  if (is.null(track) & is.null(data)) stop('must provide either a data, or both a track and a res')
  if (is.null(data) & !is.null(res) && !all(sapply(res, is.matrix))) require(raster) 
  # in case res is a list of rasters
  fixed.par.value = fixed.par
  fixed.par = !is.null(fixed.par)
  fp.comb.value = fp.comb
  fp.comb = !is.null(fp.comb)
  
  results = data.frame()
  
  if (is.null(data)) {
    #we have to make the data by ourselves based on track and res
    if (!is.null(res) && (!all(sapply(res, is.matrix)) & 
                          (class(res) != 'RasterStack' & class(res) != 'RasterLayer'))) {
      stop('res must be either a list of matrices or a RasterStack')
    } 
    res.raster = (class(res) == 'RasterStack' | class(res) == 'RasterLayer') 
    # if this is true then 'res' is a raster list, else it's a matrix list
    
    case.indices = which(track[, 'obs'] == 1)
    control.indices = which(track[, 'obs'] == 0)
    case.indices.all = which(track[, 'obs'] > 0)
    #data
    r.cases = track[case.indices, 'r']
    phi.cases = track[case.indices, 'phi']
    n.simmed = length(control.indices) / length(case.indices)
    if (res.raster) {
      # transform = !is.na(proj4string(res))
      transform = FALSE
      
      track.sp.cases = SpatialPoints(track[case.indices, c('x', 'y')],
                                     proj4string = proj4string)
      #transform points to raster coordinates
      if (transform) track.sp.cases = spTransform(track.sp.cases, CRSobj = CRS(proj4string(res)))
      res.cases = extract(res, track.sp.cases)
      
      track.sp.controls = SpatialPoints(track[control.indices, c('x', 'y')],
                                        proj4string = proj4string)
      if (transform) {
        track.sp.controls = spTransform(track.sp.controls, CRSobj = CRS(proj4string(res)))
      }
      res.controls = extract(res, track.sp.controls)
      
      track.sp.cases.all = SpatialPoints(track[case.indices.all, c('x', 'y')],
                                         proj4string = proj4string)
      if (transform) {
        track.sp.cases.all = spTransform(track.sp.cases.all, CRSobj = CRS(proj4string(res)))
      }
      res.cases.all = extract(res, track.sp.cases.all)
    } else if (!is.null(res)) {
      res.cases = data.frame(lapply(X = res, FUN = extract.matrix,
                                    x = track[case.indices, 'x'],
                                    y = track[case.indices, 'y']))
      res.controls = data.frame(lapply(X = res, FUN = extract.matrix,
                                       x = track[control.indices, 'x'],
                                       y = track[control.indices, 'y']))
      res.cases.all = data.frame(lapply(X = res, FUN = extract.matrix,
                                        x = track[case.indices.all, 'x'],
                                        y = track[case.indices.all, 'y']))
    } else {
      res.cases = 0
      res.controls = 0
      res.cases.all = 0
    }
    x.cases.all = track[case.indices.all, 'x']
    y.cases.all = track[case.indices.all, 'y']
    x.controls = track[control.indices, 'x']
    y.controls = track[control.indices, 'y']
    strata.controls = track[control.indices, 't']
  } else {
    #we just have to load the data in from the list. We need these names otherwise it won't work!
    r.cases = data$r.cases
    phi.cases = data$phi.cases
    n.simmed = data$n.simmed
    res.cases = data$res.cases
    res.controls = data$res.controls
    res.cases.all = data$res.cases.all
    x.cases.all = data$x.cases
    y.cases.all = data$y.cases
    x.controls = data$x.controls
    y.controls = data$y.controls
    strata.controls = data$strata.controls
  }
  res.cases = as.matrix(res.cases)
  res.controls = as.matrix(res.controls)
  res.cases.all = as.matrix(res.cases.all)
  res.cases[is.na(res.cases)] = 0
  res.controls[is.na(res.controls)] = 0
  res.cases.all[is.na(res.cases.all)] = 0
  #with rasters these points may become NA on occasion
  
  # ### FOR DEBUGGING PURPOSES ONLY
  # write.csv(r.cases, 'test_r_cases.csv')
  # write.csv(phi.cases, 'test_phi_cases.csv')
  # write.csv(res.cases, 'test_res_cases.csv')
  # write.csv(res.controls, 'test_res_controls.csv')
  # write.csv(res.cases.all, 'test_res_cases_all.csv')
  # write.csv(x.cases.all, 'test_x_cases_all.csv')
  # write.csv(y.cases.all, 'test_y_cases_all.csv')
  # write.csv(x.controls, 'test_x_controls.csv')
  # write.csv(y.controls, 'test_y_controls.csv')
  # write.csv(strata.controls, 'test_strata_controls.csv')
  # stop('done writing files')
  
  if (n.simmed %% 1 != 0) {
    stop('There is not an equal number of controls for every case. Check your data!')
  }
  
  if (is.null(true.values)) {
    true.values = list('null' = 0, 'resource' = 0, 'memory (1-peak)' = 0,
                       'memory (n-peak)' = 0, 'memory (power)' = 0,
                       'combination (1-peak)' = 0, 'combination (n-peak)' = 0,
                       'combination (power)' = 0)
  }
  
  #parameters
  mean.steplength = mean(r.cases)
  if (is.null(delta)) {
    #if we don't specify already
    delta = mean(r.cases < mu.stationary)
    delta = c(delta, 1-delta)
  }
  kappa.turningangle = mle.vonmises(phi.cases)$kappa
  beta.res.init = rep(0, ncol(res.cases))
  burnin = min(strata.controls) # hypothetically at what time did we start fitting the data
  
  #this code must be ran in a directory where these files exist already!
  if (model == 'null') {
    null.file = 'grizzly_null_20200528'
    # Compile the C++ file (will take some time on the first run but after that the .o and .so files are made)
    compile(paste0(null.file, '.cpp'))
    # Get the data ready
    data.null = list('r_cases' = r.cases, 'phi_cases' = phi.cases, 'delta' = delta,
                     'sl_stationary' = mu.stationary)
    # Include the parameters (with initial values - these don't mean anything because they can be changed in fit.SSF)
    par.null = list('steplength' = mean.steplength, 'angularconcentration' = kappa.turningangle,
                    'logit_lambda' = 0, 'logit_gamma' = 0)
    dyn.load(dynlib(null.file))
    # TMB conducts automatic differentiation and makes an efficient likelihood function + gradient that is called in C++
    return(MakeADFun(data = data.null, parameters = par.null, DLL = null.file))
  }
  
  if (model == 'resource') {
    res.file = 'grizzly_resource_20200616'
    compile(paste0(res.file, '.cpp'))
    
    data.res = list('n_simmed' = n.simmed, 'r_cases' = r.cases, 'phi_cases' = phi.cases,
                    'res_cases' = res.cases, 'res_controls' = res.controls, 
                    'delta' = delta, 'sl_stationary' = mu.stationary)
    par.res = list('steplength' = mean.steplength, 'angularconcentration' = kappa.turningangle,
                   'beta_res' = beta.res.init, 'logit_lambda' = 0, 'logit_gamma' = 0)
    dyn.load(dynlib(res.file))
    return(MakeADFun(data = data.res, parameters = par.res, DLL = res.file))
  }
  
  if (model == 'memory' | model == 'memory (1-peak)') {
    if (fixed.par) {
      tau1pk.file = 'grizzly_1pkmemfixed_20200624'
      compile(paste0(tau1pk.file, '.cpp'))
      data.tau1pk = list('n_simmed' = n.simmed, 'r_cases' = r.cases, 'phi_cases' = phi.cases,
                         'x_cases' = x.cases.all, 'y_cases' = y.cases.all, 
                         'x_controls' = x.controls, 'y_controls' = y.controls, 
                         'strata_controls' = strata.controls, 'delta' = delta, 
                         'sl_stationary' = mu.stationary, 'dist_coef' = fixed.par.value)
      par.tau1pk = list('steplength' = mean.steplength, 
                        'angularconcentration' = kappa.turningangle, 
                        'mean_tau' = burnin / 2, 'sd_tau' = burnin / 2, 'logit_lambda' = 0, 
                        'logit_gamma' = 0, 'log_alpha' = 0)
    } else {
      tau1pk.file = ifelse(rm.s, 'grizzly_1pkmemns_20200611', 'grizzly_1pkmem_20200805')
      compile(paste0(tau1pk.file, '.cpp'))
      data.tau1pk = list('n_simmed' = n.simmed, 'r_cases' = r.cases, 'phi_cases' = phi.cases,
                         'x_cases' = x.cases.all, 'y_cases' = y.cases.all, 
                         'x_controls' = x.controls, 'y_controls' = y.controls, 
                         'strata_controls' = strata.controls, 'delta' = delta, 
                         'sl_stationary' = mu.stationary)
      par.tau1pk = list('steplength' = mean.steplength, 
                        'angularconcentration' = kappa.turningangle, 'dist_coef' = 0, 
                        'mean_tau' = burnin / 2, 'sd_tau' = burnin / 2,
                        'logit_lambda' = 0, 'logit_gamma' = 0, 'log_alpha' = 0)
    }
    dyn.load(dynlib(tau1pk.file))
    return(MakeADFun(data = data.tau1pk, parameters = par.tau1pk, DLL = tau1pk.file))
  }
  
  if (model == 'combination' | model == 'combination (1-peak)') {
    if (fixed.par) {
      # Not used in analysis (only for testing)
      comb1pk.file = 'grizzly_1pkcombfixed_20200610'
      compile(paste0(comb1pk.file, '.cpp'))
      data.comb1pk = list('n_simmed' = n.simmed, 'r_cases' = r.cases, 'phi_cases' = phi.cases,
                          'x_cases' = x.cases.all, 'y_cases' = y.cases.all, 'x_controls' = x.controls,
                          'y_controls' = y.controls, 'strata_controls' = strata.controls, 
                          'res_casesfit' = res.cases, 'res_casesall' = res.cases.all, 
                          'res_controls' = res.controls, 'delta' = delta, 'sl_stationary' = mu.stationary, 
                          'beta.res.fixed' = fixed.par.value)
      par.comb1pk = list('steplength' = mean.steplength, 'angularconcentration' = kappa.turningangle,
                         'threshold' = 0, 'beta_res' = beta.res.init[-1], 'dist_coef' = 0, 
                         'mean_tau' = burnin / 2, 'sd_tau' = burnin / 2, 'logit_lambda' = 0, 
                         'logit_gamma' = 0, 'log_alpha' = 0)
    } else {
      comb1pk.file = 'grizzly_1pkcomb_20200717'
      compile(paste0(comb1pk.file, '.cpp'))
      
      data.comb1pk = list('n_simmed' = n.simmed, 'r_cases' = r.cases, 'phi_cases' = phi.cases,
                          'x_cases' = x.cases.all, 'y_cases' = y.cases.all, 'x_controls' = x.controls,
                          'y_controls' = y.controls, 'strata_controls' = strata.controls, 
                          'res_casesfit' = res.cases, 'res_casesall' = res.cases.all, 
                          'res_controls' = res.controls,'delta' = delta, 'sl_stationary' = mu.stationary)
      par.comb1pk = list('steplength' = mean.steplength, 'angularconcentration' = kappa.turningangle,
                         'threshold' = 0, 'beta_res' = beta.res.init, 'dist_coef' = 0, 'mean_tau' = burnin / 2,
                         'sd_tau' = burnin / 2, 'logit_lambda' = 0, 'logit_gamma' = 0, 'log_alpha' = 0)
    }
    dyn.load(dynlib(comb1pk.file))
    return(MakeADFun(data = data.comb1pk, parameters = par.comb1pk, DLL = comb1pk.file))
  }
  
  stop('Invalid "model" command. Please select from the valid commands.')
  
}

fix.raw.data = function(df.raw, proj4string = CRS('+proj=longlat +ellps=GRS80'),
                        proj4string.xy = CRS('+proj=utm +zone=8 +datum=NAD83'),
                        jitter = 0, print = TRUE) {
  #Take in raw data (will not work if this is not in the form of the grizzly bear data), and construct a time series that can be used for model fitting
  
  #df.raw - the raw data.frame with the movement tracks (time, x, y) as read in directly from CSV
  #proj4string - coordinate reference system for input data (long-lat)
  #proj4string.xy - coordinate reference system for input data (meters)
  #jitter - offset true x and y locations by runif(0, jitter). This is important because when the animal stays completely still (i.e., the same x and y coordinates) for two consecutive time steps, this produces very odd results with the memory formulation (i.e., distance from previously visited point = 0 exactly). jitter is meant to be a very small value.
  #print - give updates on the process?
  
  #return - a track that can be passed into sim.controls
  require(raster)
  
  BEARID = df.raw$BEARID[1] # will need for later
  
  if (print) print('Beginning initial data formatting')
  
  #make the dates into actual dates. Note the file needs these headings and this date format!
  df.year = as.numeric(substr(df.raw$Date, 1, 4))
  df.month = substr(df.raw$Date, 5, 6)
  df.day = substr(df.raw$Date, 7, 8)
  df.time = as.character(df.raw$Time)
  df.datetime  = as.POSIXct(paste0(df.year, '-', df.month, '-', df.day, ' ', df.time))
  df.dt = c(0, difftime(df.datetime[-1], df.datetime[-length(df.datetime)], units = 'hours'))
  df.dt = round(df.dt / 4) * 4 # round based on fix.tol
  if (any(df.dt[-1] == 0)) stop("Fix the raw data. There's a dt == 0")
  for (t in 2:length(df.datetime)) {
    #have to do this in a for loop because we need to do each one at a time
    df.datetime[t] = df.datetime[t-1] + df.dt[t] * 60 * 60 # multiply by 3600 seconds in 1 hr
  }
  
  df.code = c(1 + (df.year[-length(df.year)] < df.year[-1]), 1)
  #GUIDE FOR CODES:
  #### -1 = interpolated (bear is in the den, not moving)
  ####  0 = interpolared (moving, filling a gap in the GPS data)
  ####  1 = true, observed, with another observed point before it
  ####  2 = true, observed, with another observed point before it, last point before den entry
  ####  3 = true, observed, missing the point before it
  ####  4 = true, observed, missing the point before it, last point before den entry
  # Of course we do not want to fit interpolated points or include them in the model. This is why we generate this system. But it allows us to have a full time series without any missing data (and the row indices now indicate time since it is constant).
  
  if (print) print('Converting data to correct coordinate reference system')
  
  #convert data from degrees to meters
  df.sp = SpatialPoints(df.raw[, c('Lon', 'Lat')], proj4string = proj4string)
  df.sp.revised = spTransform(df.sp, CRSobj = proj4string.xy)@coords
  
  df.result = data.frame(t = 0:(nrow(df.raw) - 1), x = df.sp.revised[,'Lon'], 
                         y = df.sp.revised[,'Lat'], code = df.code, 
                         datetime = df.datetime, dt = df.dt)
  
  #jitter points by a very small random factor so the bear is never at the exact same point
  jitter.r = runif(nrow(df.result), 0, jitter)
  jitter.phi = runif(nrow(df.result), 0, 2*pi)
  
  df.result$x = df.result$x + jitter.r * cos(jitter.phi)
  df.result$y = df.result$y + jitter.r * sin(jitter.phi)
  
  #manually take out some of the bad data (GPS malfunctions caused some of the collars to print erroneous data after these dates)
  if (BEARID == 1008) df.result = df.result[df.result$datetime < as.POSIXct('2005-07-01 14:00:00'), ]
  if (BEARID == 1016) df.result = df.result[df.result$datetime < as.POSIXct('2005-07-01 14:00:00'), ]
  if (BEARID == 1041) df.result = df.result[df.result$datetime < as.POSIXct('2005-07-01 14:00:00'), ]
  if (BEARID == 1089) df.result = df.result[df.result$datetime < as.POSIXct('2006-07-01 14:00:00'), ]
  if (BEARID == 1107) df.result = df.result[df.result$datetime < as.POSIXct('2006-07-01 14:00:00'), ]
  if (BEARID == 1108) df.result = df.result[df.result$datetime < as.POSIXct('2006-07-01 14:00:00'), ]
  if (BEARID == 1141) df.result = df.result[df.result$datetime < as.POSIXct('2006-09-03 00:00:00'), ]
  if (BEARID == 1009) df.result = df.result[df.result$datetime < as.POSIXct('2006-01-01 00:00:00'), ]
  if (BEARID == 1043) df.result = df.result[df.result$datetime < as.POSIXct('2003-10-01 00:00:00'), ]
  if (BEARID == 1046) df.result = df.result[df.result$datetime < as.POSIXct('2006-07-01 14:00:00'), ]
  
  if (print) print('Interpolating missed GPS points')
  
  while(any(df.result$dt > 4)) {
    #interpolate points here if need be
    indices.end = which(df.result$dt > 4) # end point of interpolations
    indices.start = indices.end - 1
    codes.start = df.result$code[indices.start]
    in.den = codes.start > 0 & codes.start %% 2 == 0 #either 2 or 4
    #our interpolation for the point before indices.end is a weighted average of start & end
    weights.end = 1 - 1 / (df.result$dt[indices.end] / 4) # number of 4-hour intervals
    weights.end = weights.end * (1 - in.den) # if it's in the den we don't take the average
    #in this case weights.end = 0 because it stays in the den
    new.lat = (1 - weights.end) * df.result$y[indices.start] + weights.end * df.result$y[indices.end]
    new.lon = (1 - weights.end) * df.result$x[indices.start] + weights.end * df.result$x[indices.end]
    new.times = df.result$datetime[indices.end] - 4 * 60 * 60 # 4 hrs. before indices.end
    new.codes = (df.result$code[indices.start] != 2 & df.result$code[indices.start] != 4) - 1
    #-1 if it was in the den, 0 otherwise
    new.pts = data.frame(t = indices.end - 1 - 0.5, x = new.lon, y = new.lat, code = new.codes,
                         datetime = new.times, dt = df.result$dt[indices.end] - 4)
    df.result$dt[indices.end] = 4 # now it is!
    df.result = rbind(df.result, new.pts)
    df.result = df.result[order(df.result$t), ]
    df.result$t = 0:(nrow(df.result) - 1)
  }
  
  #now dt should be 4 everywhere
  #change some codes to 3 and 4 if necessary
  
  #where was the code 0 before and then turned to greater than 0?
  end.interpolations = df.result$code[-1] > 0 & df.result$code[-nrow(df.result)] == 0
  df.result$code[which(end.interpolations) + 1] = df.result$code[which(end.interpolations) + 1] + 2
  df.result$code[1] = 0 #we aren't fitting the first point b/c we have no step length
  
  if (print) print('Generating movement metrics for each step')
  
  #get step lengths, bearings, and turning angles
  df.steps = df.result[-1, c('x', 'y')] - df.result[-nrow(df.result), c('x', 'y')]
  df.r = c(0, apply(df.steps, 1, norm, type = '2'))
  df.theta = c(0, atan2(df.steps$y , df.steps$x))
  df.theta[is.na(df.theta)] = 0 #it didn't go anywhere so we assume the bearing is 0
  df.phi = c(0, df.theta[-1] - df.theta[-length(df.theta)]) %% (2*pi)
  
  return(data.frame(t = 0:(length(df.r)-1), x = df.result$x, y = df.result$y, 
                    theta = df.theta, phi = df.phi, r = df.r, 
                    code = df.result$code, datetime = df.result$datetime))
  
}

make.SSF.raw.data = function(track, res = NULL, write = NULL, res.bounds = NULL,
                             proj4string = as.character(NA)) {
  #track - GPS track of animal locations (returned by fix.raw.data potentially)
  #res - Raster object or list of matrices representing resources - if NULL leave resources out
  #res.bounds - lower and upper bounds (Julian date) for resource availability (2xn matrix)
  #write - file directory to place all the data files in only if we want
  #proj4string - the projection associated with the rasters
  
  # Incorporates resources into existing animal telemetry data to produce files necessary for fitting models (see generateFunction)
  
  if (is.null(data) & !is.null(res) && !all(sapply(res, is.matrix))) require(raster) 
  # in case res is a list of rasters
  
  results = data.frame()
  
  #we have to make the data by ourselves based on track and res
  if (!is.null(res) && ((class(res) != 'RasterStack' & class(res) != 'RasterLayer')) && (!all(sapply(res, is.matrix)))) {
    stop('res must be either a list of matrices or a RasterStack')
  } 
  res.raster = (class(res) == 'RasterStack' | class(res) == 'RasterLayer') # if this is true then 'res' is a raster list, else it's a matrix list
  
  case.indices = which(track[, 'obs'] == 1)
  control.indices = which(track[, 'obs'] == 0)
  case.indices.all = which(track[, 'obs'] > 0)
  #Julian dates for resource bounding (if necessary)
  jdays.cases = as.POSIXlt(track[case.indices, 'datetime'])$yday
  jdays.controls = as.POSIXlt(track[control.indices, 'datetime'])$yday
  jdays.cases.all = as.POSIXlt(track[case.indices.all, 'datetime'])$yday
  #data
  r.cases = track[case.indices, 'r']
  phi.cases = track[case.indices, 'phi']
  n.simmed = length(control.indices) / length(case.indices)
  if (res.raster) {
    transform = !is.na(proj4string(res))
    
    track.sp.cases = SpatialPoints(track[case.indices, c('x', 'y')],
                                   proj4string = CRS(proj4string))
    #transform points to raster coordinates
    if (transform) track.sp.cases = spTransform(track.sp.cases, CRSobj = CRS(proj4string(res)))
    res.cases = extract(res, track.sp.cases)
    
    track.sp.controls = SpatialPoints(track[control.indices, c('x', 'y')],
                                      proj4string = CRS(proj4string))
    if (transform) track.sp.controls = spTransform(track.sp.controls, CRSobj = CRS(proj4string(res)))
    res.controls = extract(res, track.sp.controls)
    
    track.sp.cases.all = SpatialPoints(track[case.indices.all, c('x', 'y')],
                                       proj4string = CRS(proj4string))
    if (transform) track.sp.cases.all = spTransform(track.sp.cases.all, CRSobj = CRS(proj4string(res)))
    res.cases.all = extract(res, track.sp.cases.all)
  } else if (!is.null(res)) {
    res.cases = data.frame(lapply(X = res, FUN = extract.matrix,
                                  x = track[case.indices, 'x'],
                                  y = track[case.indices, 'y']))
    res.controls = data.frame(lapply(X = res, FUN = extract.matrix,
                                     x = track[control.indices, 'x'],
                                     y = track[control.indices, 'y']))
    res.cases.all = data.frame(lapply(X = res, FUN = extract.matrix,
                                      x = track[case.indices.all, 'x'],
                                      y = track[case.indices.all, 'y']))
  } else {
    res.cases = 0
    res.controls = 0
    res.cases.all = 0
  }
  x.cases.all = track[case.indices.all, 'x']
  y.cases.all = track[case.indices.all, 'y']
  x.controls = track[control.indices, 'x']
  y.controls = track[control.indices, 'y']
  strata.controls = track[control.indices, 't']
  
  res.cases = as.matrix(res.cases)
  res.controls = as.matrix(res.controls)
  res.cases.all = as.matrix(res.cases.all)
  
  #Revise resource data for the correct temporal bounds
  if (is.null(res.bounds)) {
    res.bounds = matrix(nrow = 2, ncol = ncol(res.cases), data = c(-1,366))
  }
  for (r in 1:ncol(res.cases)) {
    res.cases[jdays.cases < res.bounds[1, r], r] = 0
    res.cases[jdays.cases > res.bounds[2, r], r] = 0
    res.controls[jdays.controls < res.bounds[1, r], r] = 0
    res.controls[jdays.controls > res.bounds[2, r], r] = 0
    res.cases.all[jdays.cases.all < res.bounds[1, r], r] = 0
    res.cases.all[jdays.cases.all > res.bounds[2, r], r] = 0
  }
  
  if (!is.null(write)) {
    #write all files to CSV in the same folder
    write.csv(r.cases, paste0(write, '/r_cases.csv'))
    write.csv(phi.cases, paste0(write, '/phi_cases.csv'))
    write.csv(n.simmed, paste0(write, '/n_simmed.csv'))
    write.csv(res.cases, paste0(write, '/res_cases.csv'))
    write.csv(res.controls, paste0(write, '/res_controls.csv'))
    write.csv(res.cases.all, paste0(write, '/res_cases_all.csv'))
    write.csv(x.cases.all, paste0(write, '/x_cases_all.csv'))
    write.csv(y.cases.all, paste0(write, '/y_cases_all.csv'))
    write.csv(x.controls, paste0(write, '/x_controls.csv'))
    write.csv(y.controls, paste0(write, '/y_controls.csv'))
    write.csv(strata.controls, paste0(write, '/strata_controls.csv'))
  }
  
  return(list(r.cases = r.cases, phi.cases = phi.cases, n.simmed = n.simmed,
              res.cases = res.cases, res.controls = res.controls, 
              res.cases.all = res.cases.all, x.cases.all = x.cases.all,
              y.cases.all = y.cases.all, x.controls = x.controls,
              y.controls = y.controls, strata.controls = strata.controls))
}

read.SSF.raw.data = function(dir) {
  #dir - file directory with all the files
  
  r.cases = read.csv(paste0(dir, '/r_cases.csv'))[, -1]
  phi.cases = read.csv(paste0(dir, '/phi_cases.csv'))[, -1]
  #r.controls = read.csv(paste0(dir, '/r_controls.csv'))[, -1]
  #phi.controls = read.csv(paste0(dir, '/phi_controls.csv'))[, -1]
  n.simmed = read.csv(paste0(dir, '/n_simmed.csv'))[, -1]
  res.cases = read.csv(paste0(dir, '/res_cases.csv'))[, -1]
  res.controls = read.csv(paste0(dir, '/res_controls.csv'))[, -1]
  res.cases.all = read.csv(paste0(dir, '/res_cases_all.csv'))[, -1]
  x.cases.all = read.csv(paste0(dir, '/x_cases_all.csv'))[, -1]
  y.cases.all = read.csv(paste0(dir, '/y_cases_all.csv'))[, -1]
  x.controls = read.csv(paste0(dir, '/x_controls.csv'))[, -1]
  y.controls = read.csv(paste0(dir, '/y_controls.csv'))[, -1]
  strata.controls = read.csv(paste0(dir, '/strata_controls.csv'))[, -1]
  
  return(list(r.cases = r.cases, phi.cases = phi.cases, n.simmed = n.simmed,
              res.cases = res.cases, res.controls = res.controls, 
              res.cases.all = res.cases.all, x.cases.all = x.cases.all,
              y.cases.all = y.cases.all, x.controls = x.controls,
              y.controls = y.controls, strata.controls = strata.controls))
  
}

write.SSF.raw.data = function(data, dir.out) {
  #Writes SSF data for exporting elsewhere to CSV. Note that make.SSF.raw.data can do this also.
  
  #data - data in the form of read.SSF.raw.data
  #dir.out - directory (folder) where files will go
  
  if (!dir.exists(dir.out)) dir.create(dir.out)
  
  write.csv(data$r.cases, paste0(dir.out, '/r_cases.csv'))
  write.csv(data$phi.cases, paste0(dir.out, '/phi_cases.csv'))
  write.csv(data$n.simmed, paste0(dir.out, '/n_simmed.csv'))
  write.csv(data$res.cases, paste0(dir.out, '/res_cases.csv'))
  write.csv(data$res.controls, paste0(dir.out, '/res_controls.csv'))
  write.csv(data$res.cases.all, paste0(dir.out, '/res_cases_all.csv'))
  write.csv(data$x.cases.all, paste0(dir.out, '/x_cases_all.csv'))
  write.csv(data$y.cases.all, paste0(dir.out, '/y_cases_all.csv'))
  write.csv(data$x.controls, paste0(dir.out, '/x_controls.csv'))
  write.csv(data$y.controls, paste0(dir.out, '/y_controls.csv'))
  write.csv(data$strata.controls, paste0(dir.out, '/strata_controls.csv'))
}

cut.raw.data = function(data, begin, end, train.begin = 0) {
  # Takes a fraction of an existing data set for model fitting (used to analyze bears year-by-year)
  
  #data - a list of raw data, in the form of something returned by read.SSF.raw.data
  #begin - in [0,1), the beginning of the sequence we want
  #end - in (begin,1], the end of the sequence we want
  #train.begin - lowest strata to include in bear's knowledge (i.e., no matter how many points are included for fitting we can also change how many points are included for "training" memory)
  
  if (end <= begin) stop("'end' must be bigger than 'begin'")
  
  n.cases = length(data$r.cases)
  n.simmed = data$n.simmed
  begin.index = ceiling(begin * n.cases)
  if (begin.index == 0) begin.index = 1
  end.index = ceiling(end * n.cases)
  if (begin.index == end.index) stop('provide a bigger interval please')
  begin.index.ctrl = (begin.index - 1) * n.simmed + 1
  end.index.ctrl = end.index * n.simmed
  
  if (data$strata.controls[begin.index.ctrl] <= train.begin) {
    stop("can't begin training after data has started")
  }
  
  train.indices = (train.begin):(max(data$strata.controls)) + 1
  new.strata = data$strata.controls[begin.index.ctrl:end.index.ctrl] - train.begin
  
  return(list(r.cases = data$r.cases[begin.index:end.index], 
              phi.cases = data$phi.cases[begin.index:end.index], 
              r.controls = data$r.controls[begin.index.ctrl:end.index.ctrl],
              phi.controls = data$phi.controls[begin.index.ctrl:end.index.ctrl],
              n.simmed = n.simmed, res.cases = data$res.cases[begin.index:end.index, ],
              res.controls = data$res.controls[begin.index.ctrl:end.index.ctrl, ], 
              res.cases.all = data$res.cases.all[train.indices, ], 
              x.cases.all = data$x.cases.all[train.indices],
              y.cases.all = data$y.cases.all[train.indices], 
              x.controls = data$x.controls[begin.index.ctrl:end.index.ctrl],
              y.controls = data$y.controls[begin.index.ctrl:end.index.ctrl],
              strata.controls = new.strata))
}

subset.controls = function(dir.in, dir.out, n.controls.new) {
  #Make a new directory with a subset of the simulated points from original data
  
  #dir.in - folder with all the raw data
  #dir.out - where we want to put all the new raw data
  #n.controls.new - number of controls to include in the real data
  
  if (!dir.exists(dir.out)) dir.create(dir.out)
  data.raw.all = read.SSF.raw.data(dir.in)
  
  n.controls.old = data.raw.all$n.simmed
  n.cases = length(data.raw.all$r.cases)
  control.keep.indices = rep(1:n.controls.new, n.cases) + 
    rep(((1:n.cases)-1)*n.controls.old, each = n.controls.new)
  
  write.csv(data.raw.all$r.cases, paste0(dir.out, '/r_cases.csv'))
  write.csv(data.raw.all$r.controls[control.keep.indices],
            paste0(dir.out, '/r_controls.csv'))
  write.csv(data.raw.all$phi.cases, paste0(dir.out, '/phi_cases.csv'))
  write.csv(data.raw.all$phi.controls[control.keep.indices],
            paste0(dir.out, '/phi_controls.csv'))
  write.csv(n.controls.new, paste0(dir.out, '/n_simmed.csv'))
  write.csv(data.raw.all$res.cases, paste0(dir.out, '/res_cases.csv'))
  write.csv(data.raw.all$res.controls[control.keep.indices, ], 
            paste0(dir.out, '/res_controls.csv'))
  write.csv(data.raw.all$res.cases.all, paste0(dir.out, '/res_cases_all.csv'))
  write.csv(data.raw.all$x.cases.all, paste0(dir.out, '/x_cases_all.csv'))
  write.csv(data.raw.all$y.cases.all, paste0(dir.out, '/y_cases_all.csv'))
  write.csv(data.raw.all$x.controls[control.keep.indices], 
            paste0(dir.out, '/x_controls.csv'))
  write.csv(data.raw.all$y.controls[control.keep.indices], 
            paste0(dir.out, '/y_controls.csv'))
  write.csv(data.raw.all$strata.controls[control.keep.indices], 
            paste0(dir.out, '/strata_controls.csv'))
}

plot.lik.surface = function(fun, gridsize, lower, upper, col = tim.colors(64), plot = TRUE, ...) {
  #A flexible function that plots a heatmap of any function mapping from R x R -> R (R is all real numbers). Can be used to plot a 2-D cross-section of the likelihood surface by making fun a customized version of a function returned from generateFunction. Used mostly for testing purposes.
  
  #fun - a function with two numeric arguments that returns a single numeric value
  #gridsize - size of the grid to produce
  #lower - 2-vector of lower bounds for plotting
  #upper - 2-vector of upper bounds for plotting
  #plot - do we make the plot (either way, returns the matrix being plotted)
  #... - additional arguments to image
  
  require(fields) #image.plot
  
  m = gridsize + 1
  xx = seq(lower[1], upper[1], length.out = m)
  yy = seq(lower[2], upper[2], length.out = m)
  
  xx.apply = rep(xx, each = m)
  yy.apply = rep(yy, times = m)
  pts = mapply(fun, xx.apply, yy.apply)
  if (plot & (any(is.na(pts)) || any(is.infinite(pts)))) stop('some points were NA or Inf')
  M = matrix(nrow = m, data = pts, byrow = TRUE)
  if (plot) {
    par(mar = c(5.1, 4.1, 4.1, 5.1))
    image(M, axes = FALSE, col = col, ...)
    axis(2, at = seq(0, 1, 0.2), seq(lower[2], upper[2], length.out = 6)) # custom y-axis
    axis(1, at = seq(0, 1, 0.2), seq(lower[1], upper[1], length.out = 6)) # custom x-axis
    image.plot(M, col = col, legend.only = TRUE, legend.mar = 4)
  }
  M
}

test.lb.sdtau.tracks = function(N, gridsize = 100, burnin = 250, n.controls = 10,
                                init.cond = 1, proj4string = CRS(as.character(NA)),
                                filename = NULL, res.layers = NULL, write.iter = FALSE,
                                lb.sdtau = 1, ...) {
  # Used to test the efficacy of different lower bounds on sigma using simulated data.
  
  #N - number of tracks to simulate and fit
  #gridsize - size of grid
  #burnin - when to start fitting points
  #n.controls - number of points to simulate
  #init.cond - number of initial parameter conditions to test in memory models
  #proj4string - projection for raster points
  #filename - if NULL, don't write the resulting data.frame to a file, otherwise the desired file name
  #res.layers - if NULL, make layers with NLMR; otherwise, read layers in
  #write.iter - should we overwrite the same file every simulation?
  #lb.sdtau - lower bound (or vector of lower bounds) for sd.tau in optimization routine
  #... - additional arguments to sim.track
  
  results.total = data.frame()
  for (sim in 1:N) {
    if (!is.null(res.layers)) {
      layers.sim = paste0('NLMs/NLM_', gridsize, '_n',(1+(sim-1)*3):(3*sim), '.csv')
    } else {
      layers.sim = NULL
    }
    trk = sim.track(gridsize = gridsize, n.res = 3, burnin = burnin, 
                    store.pars = TRUE, res.layers = layers.sim, ...)
    true.values = attr(trk, 'true')
    trk.tofit = sim.controls(track = trk, n.controls = n.controls, gridsize = gridsize, burnin = burnin)
    
    sd.true = true.values[['memory (1-peak)']][5] # the true value for sd.tau
    n.peaks = 1 + (length(true.values[['memory (n-peak)']]) - length(true.values[['memory (1-peak)']]))/3
    
    ind.sdtau = 1 # index for lb.sdtau
    continue.iterating = TRUE # did we get a good parameter estimate or must we keep going?
    while(continue.iterating & ind.sdtau <= length(lb.sdtau)) {
      if (class(B) == 'RasterStack' | class(B) == 'RasterLayer') {
        results.sim = fit.SSF(track = trk.tofit, res = B, true.values = true.values, 
                              n.peaks = n.peaks, init.cond = init.cond, proj4string = proj4string, 
                              models = 'memory', lb.sdtau = lb.sdtau[ind.sdtau])
      } else {
        results.sim = fit.SSF(track = trk.tofit, res = lapply(seq(dim(B)[3]), function(x) B[ , , x]), 
                              true.values = true.values, n.peaks = n.peaks, init.cond = init.cond, 
                              models = 'memory', lb.sdtau = lb.sdtau[ind.sdtau])
      }
      diff.sdtau = abs(lb.sdtau[ind.sdtau] - results.sim$estimate[5]) 
      # how far is the estimate from the lower bound
      continue.iterating = diff.sdtau < 1 # if the difference is less than 1 we keep going
      ind.sdtau = ind.sdtau + 1
    }
    
    results.total = rbind(results.total, cbind(N = sim, results.sim, lb.sdtau = lb.sdtau[ind.sdtau - 1]))
    if (write.iter & !is.null(filename)) write.csv(x = results.total, file = filename)
  }
  
  if (!is.null(filename)) write.csv(x = results.total, file = filename)
  results.total
  
}

sim.and.fit.tracks = function(N, gridsize = 100, burnin = 2, n.controls = 10, n.res = 3,
                              init.cond = 1, proj4string = CRS(as.character(NA)), 
                              models = c('null', 'resource', 'memory', 'combination'),
                              filename = NULL, res.layers = NULL, write.iter = FALSE, 
                              BIC = FALSE, same.layers = FALSE, raster = FALSE, rest = 0, 
                              mu.stationary = 0.01, parallel = 0, ...) {
  # Perform the simulation analysis from Section 3.3 in the paper, all at once. Simulates multiple tracks (N of them to be exact) and fits the prescribed models. Takes a long time - this function is best used on a server or some other form of remote computer!
  
  #N - number of tracks to simulate and fit
  #gridsize - size of grid
  #burnin - when to start fitting points
  #n.controls - number of points to simulate
  #n.res - as in sim.track
  #init.cond - number of initial parameter conditions to test in memory models
  #proj4string - projection for raster points
  #models - list of models we intend to fit here (same notation as in generateFunction)
  #filename - if NULL, don't write the resulting data.frame to a file, otherwise the desired file name
  #res.layers - if NULL, make layers with NLMR; otherwise, read layers in
  #write.iter - should we overwrite the same file every simulation?
  #BIC - use BIC or AIC?
  #same.layers - use the same resource layers for each time?
  #raster - are we using rasters?
  #rest - how long of a rest period should we include? default no rest; business as usual
  #mu.stationary - stationary step length parameter
  #parallel - number of cores to use
  #... - additional arguments to sim.track
  
  #return - a data frame of all the results for each simulated track
  
  ext = ifelse(raster, '.tif', '.csv')
  
  if (parallel == 0) {
    # No parallel - we use a for loop...
    results.total = data.frame()
    
    for (sim in 1:N) {
      if (!is.null(res.layers)) {
        if (same.layers) {
          layers.sim = paste0(res.layers, '/NLM_', gridsize, '_n', 
                              1:n.res, ext)
        } else {
          layers.sim = paste0(res.layers, '/NLM_', gridsize, '_n', 
                              (1+(sim-1)*n.res):(n.res*sim),  ext)
        }
      } else {
        layers.sim = NULL
      }
      # Simulate the track 
      trk = sim.track(gridsize = gridsize, n.res = n.res, burnin = burnin, 
                      store.pars = TRUE, res.layers = layers.sim, raster = raster,
                      mu.s = mu.stationary, ...)
      # Add resting point if necessary (it usually is not)
      trk = add.resting.point(trk, round(runif(1, burnin, nrow(trk))), rest)
      true.values = attr(trk, 'true')
      if (is.null(gridsize) & !raster) B = 0
      if (class(B) == 'RasterStack' | class(B) == 'RasterLayer') {
        # Fit the models to simulated tracks
        trk.tofit = sim.controls(track = trk, n.controls = n.controls, gridsize = gridsize, burnin = burnin, 
                                 raster = sum(B))
        results.sim = fit.SSF(track = trk.tofit, res = B, true.values = true.values, 
                              n.peaks = 1 + (length(true.values[['memory (n-peak)']]) - length(true.values[['memory (1-peak)']]))/3,
                              power = FALSE, init.cond = init.cond, proj4string = proj4string, models = models,
                              BIC = BIC, lb.sdtau = 20, mu.stationary = mu.stationary, 
                              control = list(eval.max = 500, iter.max = 500))
      } else {
        if (is.null(gridsize)) {
          res = NULL
        } else {
          res = lapply(seq(dim(B)[3]), function(x) B[ , , x])
        }
        trk.tofit = sim.controls(track = trk, n.controls = n.controls, gridsize = gridsize, burnin = burnin)
        # Fit the models to simulated tracks
        results.sim = fit.SSF(track = trk.tofit, res = res, true.values = true.values, 
                              n.peaks = 1 + (length(true.values[['memory (n-peak)']]) - length(true.values[['memory (1-peak)']]))/3,
                              power = FALSE, init.cond = init.cond, models = models, 
                              BIC = BIC, lb.sdtau = 20, mu.stationary = mu.stationary, 
                              control = list(eval.max = 500, iter.max = 500))
      }
      results.total = rbind(results.total, cbind(N = sim, results.sim))
      if (write.iter & !is.null(filename)) write.csv(x = results.total, file = filename)
    }
  } else {
    # Use parallel programming here to do each track on its own core
    require(doParallel)
    registerDoParallel(cores = parallel)
    if (write.iter) warning("'write.iter will not do anything if parallel > 0")
    results.total = foreach(sim = 1:N, .combine = rbind, .errorhandling = 'stop') %dopar% {
      if (!is.null(res.layers)) {
        if (same.layers) {
          layers.sim = paste0(res.layers, '/NLM_', gridsize, '_n', 
                              1:n.res, ext)
        } else {
          layers.sim = paste0(res.layers, '/NLM_', gridsize, '_n', 
                              (1+(sim-1)*n.res):(n.res*sim), ext)
        }
      } else {
        layers.sim = NULL
      }
      trk = sim.track(gridsize = gridsize, n.res = n.res, burnin = burnin, 
                      store.pars = TRUE, res.layers = layers.sim, raster = raster,
                      mu.s = mu.stationary, ...)
      trk = add.resting.point(trk, round(runif(1, burnin, nrow(trk))), rest)
      true.values = attr(trk, 'true')
      if (is.null(gridsize) & !raster) B = 0
      if (class(B) == 'RasterStack' | class(B) == 'RasterLayer') {
        trk.tofit = sim.controls(track = trk, n.controls = n.controls, gridsize = gridsize, burnin = burnin,
                                 raster = sum(B))
        results.sim = fit.SSF(track = trk.tofit, res = B, true.values = true.values, 
                              n.peaks = 1 + (length(true.values[['memory (n-peak)']]) - length(true.values[['memory (1-peak)']]))/3,
                              power = FALSE, init.cond = init.cond, proj4string = proj4string, models = models,
                              BIC = BIC, lb.sdtau = 20, mu.stationary = mu.stationary, 
                              control = list(eval.max = 500, iter.max = 500))
      } else {
        if (is.null(gridsize)) {
          res = NULL
        } else {
          res = lapply(seq(dim(B)[3]), function(x) B[ , , x])
        }
        trk.tofit = sim.controls(track = trk, n.controls = n.controls, gridsize = gridsize, burnin = burnin)
        results.sim = fit.SSF(track = trk.tofit, res = res, true.values = true.values, 
                              n.peaks = 1 + (length(true.values[['memory (n-peak)']]) - length(true.values[['memory (1-peak)']]))/3,
                              power = FALSE, init.cond = init.cond, models = models, 
                              BIC = BIC, lb.sdtau = 20, mu.stationary = mu.stationary, 
                              control = list(eval.max = 500, iter.max = 500))
      }
      cbind(N = sim, results.sim)
    }
  }
  
  if (!is.null(filename)) write.csv(x = results.total, file = filename)
  results.total
}

write.NLM = function(gridsize, n) {
  # Write a simulated landscape to CSV (for later use)
  
  #gridsize - integer (vector) of grid sizes
  #n - how many NLMs to write
  
  require(NLMR)
  
  for (gs in gridsize) {
    for (i in 1:n) {
      M = matrix(nrow = gs, data = nlm_gaussianfield(gs, gs)@data@values)[,gs:1]
      write.csv(M, paste0('NLMs/NLM_', gs, '_n', i, '.csv'))
    }
  }
  
  0 # return if everything went OK
}

add.par.names = function(dir) {
  # Helper function to make output data from fit.SSF look nicer. Sets the row names for the data.frame to the associated parameter using the ordering below.
  files = list.files(dir)
  
  PAR_VARIABLE_NAMES = c('mu', 'kappa', 'lambda', 'gamma', 'mu', 'kappa', 
                         'beta_1', 'beta_2', 'beta_3', 'beta_4', 'beta_5',
                         'beta_6', 'lambda', 'gamma', 'mu', 'kappa', 'beta_mem',
                         'mu_mem', 'sigma_mem', 'lambda', 'gamma', 'alpha',
                         'mu', 'kappa', 'beta_0', 'beta_1', 'beta_2', 'beta_3',
                         'beta_4', 'beta_5', 'beta_6', 'beta_mem', 'mu_mem', 
                         'sigma_mem', 'lambda', 'gamma', 'alpha')
  
  for (file in files) {
    read.file = read.csv(paste0(dir, '/', file))
    
    n.simulations = nrow(read.file) / length(PAR_VARIABLE_NAMES)
    colnames(read.file)[1] = 'Par'
    read.file$Par = rep(PAR_VARIABLE_NAMES, n.simulations)
    
    # colnames(read.file)[2] = 'N'
    # read.file$N = rep(1:n.simulations, each = length(PAR_VARIABLE_NAMES))
    
    write.csv(read.file, paste0(dir, '/', file)) # overwrite with same name
  }
  
  0 # return if everything went OK
}

par.bias.MSE = function(file, mod) {
  #Return bias and MSE measurements for a given parameter
  
  #file - results data.frame filename
  #mod - model (long names in files) to get parameters for
  
  result = data.frame()
  
  par.estimates = read.csv(file)
  
  this.model = subset(par.estimates, model == mod)
  n.simulations = length(unique(par.estimates$N))
  n.pars = nrow(this.model) / n.simulations
  par.names = unique(this.model$Par)
  
  for (p in 1:n.pars) {
    this.par = subset(this.model, Par == par.names[p])
    true.val = this.par[1,'true']
    
    if (par.names[p] == 'alpha' | par.names[p] == 'lambda' | 
        par.names[p] == 'beta_0' | par.names[p] == 'beta_d') {
      #transform back
      this.par$estimate = 1/(1+exp(-this.par$estimate))
      true.val = 1/(1+exp(-true.val))
    }
    if (par.names[p] == 'gamma') {
      this.par$estimate = this.par$estimate / log(10)
      true.val = true.val / log(10)
    }
    
    result = rbind(result, c(par.name = par.names[p], true.value = true.val, 
                             bias = mean(this.par$estimate) - true.val,
                             MSE = mean((this.par$estimate - true.val)^2), 
                             N = n.simulations))
  }
  
  result
}

change.par.names = function(dir, n.res,
                            mods = c('null', 'resource', 'memory', 'combination')) {
  #Changes the first column to the proper parameter names (the quick and dirty way)
  
  #dir - file directory with results data.frames
  #n.res - number of resource parameters
  #mods - list of models
  
  par.names = c()
  
  for (mod in mods) {
    if (mod == 'null') {
      par.names = c(par.names, 'rho_ns', 'kappa', 'lambda', 'gamma')
    } else if (mod == 'resource') {
      par.names = c(par.names, 'rho_ns', 'kappa', paste0('beta', 1:n.res),
                    'lambda', 'gamma')
    } else if (mod == 'memory') {
      par.names = c(par.names, 'rho_ns', 'kappa', 'beta_d', 'mu', 'sigma',
                    'alpha', 'lambda', 'gamma')
    } else if (mod == 'combination') {
      par.names = c(par.names, 'rho_ns', 'kappa', 'beta_0',
                    paste0('beta', 1:n.res), 'beta_d', 'mu', 'sigma',
                    'alpha', 'lambda', 'gamma')
    }
  }
  
  n.per.sim = length(par.names)
  
  files = list.files(dir)
  for (filename in files) {
    file = read.csv(paste0(dir, '/', filename))
    n.sims = nrow(file) / n.per.sim
    file[,1] = rep(par.names, n.sims)
    colnames(file)[1] = 'Par'
    write.csv(file, paste0(dir, '/', filename))
  }
  
  return(0) # if everything went OK
  
}

plot.par.estimates = function(dir, mod, par = NULL, ylim.tol = 5000, type = c('each', 'sep', 'hist'), ...) {
  #dir - file directory with results data.frames
  #mod - exact model we want to get parameters for
  #par - specific parameter to graph
  #ylim.tol - maximum confidence interval width to 
  #type - do we do a histogram or what?
  #... - additional arguments to boxplot
  
  #returns - a boxplot of parameter estimates for each variable
  
  if (!is.null(par)) warning('could be missing parameters and giving really odd results FYI')
  if (length(type) > 1) type = type[1]
  if (type == 'sep' & is.null(par)) stop('if type == "sep" then must have "par"')
  
  hist = type == 'hist'
  
  files = list.files(dir)
  first = TRUE
  
  for(f in files) {
    results = read.csv(paste0(dir, '/', f))
    results = cbind(results, CL = results$estimate - 2*results$se, CU = results$estimate + 2*results$se)
    #I need to remind myself why this is here
    
    this.model = subset(results, model == mod)
    n.simulations = length(unique(results$N))
    
    if (type == 'sep' & first) {
      all.estimates = as.data.frame(matrix(ncol = 2, nrow = 0))
      # all.estimates = matrix(ncol = 0, nrow = n.simulations)
      first = FALSE
    }
    
    n.pars = nrow(this.model) / n.simulations # number of parameters for this model
    par.names = unique(this.model$Par)
    if (!is.null(par)) par.names = par
    
    for (p in 1:length(par.names)) {
      this.par = subset(this.model, Par == par.names[p])
      true.val = this.par[1,'TRUE.']
      if (hist) {
        this.par = this.par[, 'estimate']
        hist(this.par, main = paste0(f, ': ', par.names[p]))
        abline(v = 0, col = 'black')
        abline(v = true.val, col = 'red')
      } else {
        this.par = this.par[, c('estimate', 'CL', 'CU')]
        if (type == 'each') {
          this.par = subset(this.par, !is.na(CL) & (CU - CL < ylim.tol))
          boxplot(t(this.par), main = paste0(f, ': ', par.names[p], '\n',
                                             n.simulations - nrow(this.par), ' left out'), ...)
          abline(h = 0, col = 'black')
          abline(h = true.val, col = 'red')
        } else {
          all.estimates = rbind(all.estimates, data.frame(file = f, 
                                                          value = this.par$estimate))
          # all.estimates = cbind(all.estimates, this.par$estimate)
        }
      }
    }
  }
  
  if (type == 'sep') {
    # Transformation here depends on which parameter is being plotted
    require(ggplot2)
    # ggplot(all.estimates, aes(x=file, y=value)) + geom_violin() + 
    #   scale_x_discrete(labels = c('T = 600; K = 10', 'T = 600; K = 50',
    #                               'T = 1200; K = 10', 'T = 1200; K = 50')) + 
    #   xlab('') + ylab('Parameter estimate value') + ylim(c(-5, 1)) + 
    #   geom_hline(yintercept = log10(1/60), color = 'red')
    
    # for alpha
    all.estimates$new_value = all.estimates$value / log(10)
    ggplot(all.estimates, aes(x=file, y=new_value)) + geom_violin() +
      scale_x_discrete(labels = c('T = 600; K = 10', 'T = 600; K = 50',
                                  'T = 1200; K = 10', 'T = 1200; K = 50')) +
      xlab('') + ylab('Parameter estimate value') + ylim(c(-6, 2)) +
      geom_hline(yintercept = log10(1/60), color = 'red')
    
    # for beta_1
    # ggplot(all.estimates, aes(x=file, y=value)) + geom_violin() + 
    #   scale_x_discrete(labels = c('T = 600; K = 10', 'T = 600; K = 50',
    #                               'T = 1200; K = 10', 'T = 1200; K = 50')) + 
    #   xlab('') + ylab('Parameter estimate value') + ylim(c(0, 15)) + 
    #   geom_hline(yintercept = 7.5, color = 'red')
  }
}

AIC.table = function(dir, weights = FALSE) {
  #get an AIC table for a set of simulations
  
  #dir - filename
  #weights - do we want IC weights instead of vote count?
  
  #return a list of which model was best for each run
  if (weights) require(qpcR)
  res = list()
  
  files = list.files(dir)
  for (f in files) {
    results = read.csv(paste0(dir, '/', f))
    
    n.models = length(unique(results$model))
    n.simulations = length(unique(results$N))
    AIC.this = numeric(n.models)
    if (nrow(subset(results, N == 1)) == 31) {
      #n.res = 3
      IC.indices = c(1,5,12,31)
    } else {
      #n.res = 6
      IC.indices = c(1,5,15,37)
    }
    for (n in 1:n.simulations) {
      this.result = subset(results, N == n)
      this.AIC = this.result$IC[IC.indices]
      if (weights) {
        AIC.this = AIC.this + akaike.weights(this.AIC)$weights / n.simulations
      } else {
        AIC.ind = which(this.AIC == min(this.AIC))[1] #if there's a tie take the first one
        AIC.this[AIC.ind] = AIC.this[AIC.ind] + 1
      }
    }
    res[[f]] = AIC.this
  }
  
  res
}

getAvgDistRatio = function(mu, kappa, t, n.persim, n.sim, ...) {
  # Helper function used mostly for testing - evaluates the average distance that the simulated animal is from a given point on its track at any time
  
  #mu - mean step length
  #kappa - angular concentration
  #t - get average distance from your current location to your location t points away
  #n.persim - sample size reuqired for this average
  #n.sim - number of tracks to simulate
  #... - additional parameters to sim.track (i.e., if we want to include memory - no resources!)
  
  ratios = numeric(n.sim)
  
  for (i in 1:n.sim) {
    trk = sim.track(mu = mu, kappa = kappa, gridsize = NULL, n.steps = t + n.persim, ...)
    end.indices = (t+1):(t+n.persim+1)
    start.indices = end.indices - t
    compare.indices = start.indices + 1
    dists.curr = sqrt((trk[end.indices, 'x'] - trk[start.indices, 'x'])^2 + (trk[end.indices, 'y'] - trk[start.indices, 'y'])^2) 
    dists.comp = sqrt((trk[end.indices, 'x'] - trk[compare.indices, 'x'])^2 + (trk[end.indices, 'y'] - trk[compare.indices, 'y'])^2)
    ratios[i] = mean(abs(dists.curr - dists.comp)) / mean(dists.curr)
    print(ratios[i])
  }
  
  ratios
}

bhatt.coef = function(x, y, bw=bw.nrd0, ...) {
  # Approximates the Bhattacharyya coefficient for two numerics x and y (Not my code)
  
  #x and y are numerics; bw is the way we decide to make the bins
  #... additional arguments to bw
  if(!is.numeric(x) | !is.numeric(y)) stop('x and y must both be numeric')
  if(length(x) < 2 | length(y) < 2) stop('x and y must both have at least 2 values')
  if(length(bw) != 1) stop("'bw' must be either a single numeric value or a single function.")   
  if (!is.function(bw) & !is.numeric(bw)) stop("'bw' must be either a single numeric value or a single function.")
  if(is.numeric(bw)) bw<-round(bw)
  
  #Setting the right number of bins
  if(is.function(bw)) {
    band.width = bw(c(x,y), ...) #bin width
    #bin breaks
    minxy = min(c(x,y))
    maxxy = max(c(x,y)+band.width)
    minby = (maxxy - minxy) / 100000
    band.width = max(band.width, minby)
    bin.breaks = seq(from=minxy, to=maxxy, by=band.width) #adding an extra bandwith to the max to be sure to include all the data
    bin.n = length(bin.breaks)-1 #number of bins
  } else {
    bin.breaks = hist(c(x, y), breaks = bw, plot = FALSE)$breaks
    band.width = diff(bin.breaks)[1]
    bin.n = bw
  }
  
  #Counting the number of elements per bin
  histx = hist(x, breaks = bin.breaks, plot = FALSE)[[2]]
  histy = hist(y, breaks = bin.breaks, plot = FALSE)[[2]]
  #Relative counts
  rel.histx = histx / sum(histx)
  rel.histy = histy / sum(histy)
  
  sum(sqrt(rel.histx*rel.histy))
}

bhatt.coef.2d = function(x, y, bw=bw.nrd0, matrix = FALSE, ...) {
  require(gplots)
  
  # 2-D Bhattacharyya coefficient approximation
  
  #x and y - matrices or df's with 2 columns
  #bw - how we decide to make the bins - either a function or simply a 2-vector of breaks
  #matrix - are the inputs spatially distributed matrices or not?
  #... additional arguments to bw
  
  if (!matrix) {
    if(length(x) < 2 | length(y) < 2) stop('x and y must both have at least 2 values')
    if (!is.function(bw) & !is.numeric(bw)) stop("'bw' must be either a single numeric value or a single function.")
    if(is.numeric(bw)) bw<-round(bw)
    
    #Setting the right number of bins
    if(is.function(bw)) {
      x.all = c(x[,1],y[,1])
      band.width.1 = bw(x.all, ...) #bin width
      #bin breaks
      minx = min(x.all)
      maxx = max(x.all+band.width.1)
      #get minimum bin break in case bw returns a value too small
      minbw.1 = (maxx - minx) / 100000
      band.width.1 = max(band.width.1, minbw.1)
      bin.breaks.1 = seq(from=minx, to=maxx, by=band.width.1) #adding an extra bandwith to the max to be sure to include all the data
      
      y.all = c(x[,2],y[,2])
      band.width.2 = bw(y.all, ...) #bin width
      #bin breaks
      miny = min(y.all)
      maxy = max(y.all+band.width.2)
      #get minimum bin break in case bw returns a value too small
      minbw.2 = (maxy - miny) / 100000
      band.width.2 = max(band.width.2, minbw.2)
      bin.breaks.2 = seq(from=miny, to=maxy, by=band.width.2) #adding an extra bandwith to the max to be sure to include all the data
      
      bin.n = c(length(bin.breaks.1)-1,length(bin.breaks.2)-1)  #number of bins
    } else {
      bin.n = bw
    }
    #Counting the number of elements per bin
    histx = hist2d(x = x[,1], y = x[,2], nbins = bin.n)$counts
    histy = hist2d(x = y[,1], y = y[,2], nbins = bin.n)$counts
    
  } else {
    if (nrow(x) != nrow(y) | ncol(x) != ncol(y)) stop('x and y must have same dimensions')
    histx = x
    histy = y
    
  }
  
  #Relative counts
  rel.histx = as.numeric(histx) / sum(histx)
  rel.histy = as.numeric(histy) / sum(histy)
  
  sum(sqrt(rel.histx*rel.histy))
}

downscale.matrix = function(M, rate = 2) {
  #M - matrix
  #r - how strongly to scale
  
  #returns a matrix with dim(M)/r cells based on the average of those cells in M
  
  if (nrow(M) %% rate != 0 | ncol(M) %% rate != 0) stop('invalid rate.')
  
  row.ind = seq(1, nrow(M), rate)
  col.ind = seq(1, ncol(M), rate)
  
  row.ind = rep(row.ind, ncol(M)/rate)
  col.ind = rep(col.ind, each = nrow(M)/rate)
  
  new.data = sapply(X = 1:length(row.ind), FUN = function(x) {
    matrix.sub = mean(M[row.ind[x]:(row.ind[x]+rate-1), col.ind[x]:(col.ind[x]+rate-1)])
  })
  
  matrix(nrow = nrow(M)/rate, data = new.data)
}

prof.lik.CI = function(fun, index, optimum, gr = NULL, init.step = 1, alpha = 0.05,
                       constraints = matrix(ncol = 2, byrow = TRUE,
                                            data = rep(c(-Inf, Inf), length(optimum))), 
                       tol = 1e-3, max.steps = 200, plot = FALSE, write = NULL, 
                       custom.message = '', ...) {
  # Compute confidence intervals using likelihood profiling (see Fischer et al., 2020)
  
  #fun - negative log-likelihood function
  #index - index (1-based) of desired parameter
  #optimum - optimal value of parameters
  #gr - gradient
  #init.step - the first step length to try for the binary search
  #alpha - confidence level
  #constraints - lower and upper (respectively) optimization bounds for all parameters
  #tol - how close do we have to be to the target
  #max.steps - how long do we let the algorithm run for
  #plot - plot profile in parameter space?
  #write - write plot.pts to file?
  #custom.message - extra text to have if a new optimum is found
  #... - additional arguments to optimizer
  
  #returns - c(lower bound, upper bound)
  
  optimum.val = fun(optimum)
  target = optimum.val + qchisq(1-alpha, 1) / 2
  
  # start by going lower
  new.par = optimum[index]
  steps.lower = 0
  
  step.try = init.step
  mult.fac = 2
  
  plot.pts = cbind(value = optimum.val, t(optimum))
  
  while(abs(step.try) > tol/2 & steps.lower < max.steps) {
    while (new.par - step.try <= constraints[index, 1]) {
      if (new.par == constraints[index, 1]) {
        step.try = 0
        break
      } else {
        step.try = step.try / 2
      }
    }
    
    if (step.try != 0) {
      new.par = new.par - step.try
      pars.try = optimum
      pars.try[index] = new.par
      new.constraints = constraints
      new.constraints[index, ] = c(new.par, new.par)
      fun.opt = nlminb(pars.try, fun, gr, lower = new.constraints[, 1],
                       upper = new.constraints[, 2], ...)
      fun.val = fun.opt$objective
      plot.pts = rbind(plot.pts, c(fun.val, fun.opt$par))
      
      if (fun.val > target) {
        step.try = abs(step.try) * -0.5
        mult.fac = 0.5 # we don't need to double all the time anymore
      } else {
        # try a step of twice the size
        step.try = abs(step.try) * mult.fac
      }
      
      if (fun.val < optimum.val - 1e3) {
        warning(paste0(custom.message, 'found new optimum at: ', paste(new.par, collapse = ', ')))
      }
      
      steps.lower = steps.lower + 1
    } else {
      steps.lower = max.steps
    }
    CI.lower = new.par
  }
  
  #then go upper
  
  new.par = optimum[index]
  steps.upper = 0
  
  step.try = init.step
  mult.fac = 2
  
  while(abs(step.try) > tol/2 & steps.upper < max.steps) {
    while (new.par + step.try >= constraints[index, 2]) {
      if (new.par == constraints[index, 2]) {
        step.try = 0
        break
      } else {
        step.try = step.try / 2
      }
    }
    
    if (step.try != 0) {
      new.par = new.par + step.try
      pars.try = optimum
      pars.try[index] = new.par
      new.constraints = constraints
      new.constraints[index, ] = c(new.par, new.par)
      fun.opt = nlminb(pars.try, fun, gr, lower = new.constraints[, 1],
                       upper = new.constraints[, 2], ...)
      fun.val = fun.opt$objective
      plot.pts = rbind(plot.pts, c(fun.val, fun.opt$par))
      
      if (fun.val > target) {
        step.try = abs(step.try) * -0.5
        mult.fac = 0.5 # we don't need to double all the time anymore
      } else {
        # try a step of twice the size
        step.try = abs(step.try) * mult.fac
      }
      
      if (fun.val < optimum.val - 1e3) {
        warning(paste0(custom.message, 'found new optimum at: ', paste(new.par, collapse = ', ')))
      }
      
      steps.upper = steps.upper + 1
    } else {
      steps.upper = max.steps
    }
    CI.upper = new.par
  }
  
  if(plot) plot(plot.pts$par, plot.pts$value, xlab = 'parameter', 
                ylab = 'negative log-likelihood')
  if (!is.null(write)) write.csv(plot.pts, write)
  
  data.frame(lower = CI.lower, upper = CI.upper)
}

get.CI.lower.fromfile = function(par.index, ID, model, dir = '') {
  #ID = bear ID (for file name)
  #par.index = number of the parameter (1-based; for file name)
  #model = model in question (for file name)
  #dir = where the files are
  
  if (!model %in% c('null', 'res', 'mem', 'comb')) {
    stop('incorrect model name')
  }
  
  if (dir != '') dir = paste0(dir, '/')
  CI.file = read.csv(paste0(dir, ID, model, par.index, '.csv'))[, -1]
  
  estimate.par.value = CI.file[1, par.index + 1]
  lower.indices = which(CI.file[, par.index + 1] < estimate.par.value)
  
  if (length(lower.indices) == 0) {
    return(CI.file[1, par.index + 1])
  }
  
  return(CI.file[lower.indices[length(lower.indices)], par.index + 1])
}

get.CI.upper.fromfile = function(par.index, ID, model, dir = '') {
  #ID = bear ID (for file name)
  #par.index = number of the parameter (1-based; for file name)
  #model = model in question (for file name)
  #dir = where the files are
  
  if (!model %in% c('null', 'res', 'mem', 'comb')) {
    stop('incorrect model name')
  }
  
  if (dir != '') dir = paste0(dir, '/')
  CI.file = read.csv(paste0(dir, ID, model, par.index, '.csv'))[, -1]
  
  estimate.par.value = CI.file[1, par.index + 1]
  upper.indices = which(CI.file[, par.index + 1] > estimate.par.value)
  
  if (length(upper.indices) == 0) {
    return(CI.file[1, par.index + 1])
  }
  
  return(CI.file[upper.indices[length(upper.indices)], par.index + 1])
}

get.CI.estimate.fromfile = function(par.index, ID, model, dir = '') {
  #ID = bear ID (for file name)
  #par.index = number of the parameter (1-based; for file name)
  #model = model in question (for file name)
  #dir = where the files are
  
  if (!model %in% c('null', 'res', 'mem', 'comb')) {
    stop('incorrect model name')
  }
  
  if (dir != '') dir = paste0(dir, '/')
  CI.file = read.csv(paste0(dir, ID, model, par.index, '.csv'))[, -1]
  
  CI.file[1, par.index + 1]
}

convert.IC = function(dir, dir.out, data.size, BIC.in = TRUE) {
  #Converts BIC in results files to AIC, or the other way
  
  #dir - directory with results data frames
  #dir.out - where to write new files?
  #data.size - size of the data (for BIC calculation)
  #BIC.in - are the results currently in as BIC?
  
  files = list.files(dir)[grepl('.csv', list.files(dir))]
  if (length(data.size) == 1) data.size = rep(data.size, length(files))
  
  for (f in files) {
    file = read.csv(paste0(dir, '/', f))
    
    n.simulations = length(unique(file$N))
    n.data = data.size[which(files == f)]
    ICs = file$IC
    
    #get the total number of models
    model.list = unique(file$model[file$N == min(file$N)])
    #a vector of model parameter numbers for each model
    model.par.nums = as.numeric(table(file$model[file$N == min(file$N)])[model.list])
    #the number of times each model appears is the number of parameters
    model.par.nums = rep(rep(model.par.nums, model.par.nums), n.simulations)
    #now we have a vector the length of the data.frame
    
    if (BIC.in) {
      model.neg.log.liks = file$IC - log(n.data) * model.par.nums
      model.ICs.new = model.neg.log.liks + 2 * model.par.nums
    } else {
      model.neg.log.liks = file$IC - 2 * model.par.nums
      model.ICs.new = model.neg.log.liks + log(n.data) * model.par.nums
    }
    
    file$IC = model.ICs.new
    write.csv(file, paste0(dir.out, '/', f))
  }
  
}

get.LRT.p.values = function(file, data.size, n.simulations = NULL, BIC = TRUE) {
  #Gets a matrix of p-values for likelihood ratio tests
  
  #file - CSV file of fit results (returned from fit.SSF)
  #data.size - either a number or vector of data sizes
  #n.simulations - number of fits in the file
  #BIC - BIC or AIC?

  fit.data = read.csv(file)
  
  if (is.null(n.simulations)) n.simulations = length(unique(fit.data$N))
  if (length(data.size) == 1) data.size = rep(data.size, n.simulations)
  if (length(data.size) != n.simulations) stop('incorrect data.size length')
  
  n.rows.per.sim = nrow(fit.data) / n.simulations
  if (n.rows.per.sim %% 1 != 0) {
    stop('n.simulations may be wrong!')
  }
  
  p.values = matrix(nrow = n.simulations, ncol = 5, 1)
  # columns are for null-res, null-mem, null-comb, res-comb, mem-comb
  
  for (i in 1:n.simulations) {
    n.data = data.size[i]
    fit.data.n = fit.data[(i-1)*n.rows.per.sim + 1:n.rows.per.sim, ]
    
    model.IC.indices = sapply(X = unique(fit.data.n$model), FUN = function(x) {
      which(fit.data.n$model == x)[1]
    }) # get one index of every model fit result to get BIC
    
    model.ICs = fit.data.n$IC[model.IC.indices]
    model.par.nums = as.numeric(table(fit.data.n$model)[unique(fit.data$model)])
    #the number of times each model appears is the number of parameters
    
    if (BIC) {
      model.neg.log.liks = (model.ICs - log(n.data) * model.par.nums) / 2
    } else {
      model.neg.log.liks = (model.ICs - 2 * model.par.nums) / 2
    }
    #get actual negative log-likelihood values
    
    if (length(model.ICs) == 4) {
      #if we have 4 models
      p.null.res = 1 - pchisq(2 * (model.neg.log.liks[1] - model.neg.log.liks[2]),
                              model.par.nums[2] - model.par.nums[1])
      p.null.mem = 1 - pchisq(2 * (model.neg.log.liks[1] - model.neg.log.liks[3]),
                              model.par.nums[3] - model.par.nums[1])
      p.null.comb = 1 - pchisq(2 * (model.neg.log.liks[1] - model.neg.log.liks[4]),
                               model.par.nums[4] - model.par.nums[1])
      p.res.comb = 1 - pchisq(2 * (model.neg.log.liks[2] - model.neg.log.liks[4]),
                              model.par.nums[4] - model.par.nums[2])
      p.mem.comb = 1 - pchisq(2 * (model.neg.log.liks[3] - model.neg.log.liks[4]),
                              model.par.nums[4] - model.par.nums[3])
      
      p.values[i, ] = c(p.null.res, p.null.mem, p.null.comb, p.res.comb, p.mem.comb)
      
    } else {
      warning('What models are you fitting? Doing nothing...')
    }
    
  }
  
  p.values
}