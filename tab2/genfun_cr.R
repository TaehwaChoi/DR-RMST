## Generate a dataset that satisfies a specified censoring rate

t = Time ; p = 0.3
cens_rate = function(t, p) {
  # non-parametric estimation of survival time distribution
  dens = density(t)  # density() : returns the density data
  x = dens$x
  y = dens$y
  y.loess = loess(y ~ x, span = 0.1)
  
  # Integration 
  # censoring rate = 1-mean(delta); mean(delta) = P(T<=C)
  f_den = function(x) predict(y.loess, newdata = x)
  f_int = function(u, theta) (u/theta) * (f_den(u))
  int = function(x, theta, p) { 
        cens = integrate(f_int, theta = x, lower = min(t), upper = max(t))
        return(cens$value - p)
        }
  sol = uniroot(int, p = p, interval = c(0.001, 100))$root
  return(sol)
}

cens_rate(t = Time, p = 0.3)


