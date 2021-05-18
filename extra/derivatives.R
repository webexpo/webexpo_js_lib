source("C:\\Users\\blitzrcr\\patrickb.stat\\Lavoue\\R\\normal_distrn.R")
 
                                                    
# Version 0.2 (May 2021)
# [distributed]

# Change Log
# =================
#
# Version 0.2 (May 2021)
# ----------------------
#   Corrected an error sign in fcts
#     - d2logPhi.dmu2 
#     - d2logPhi.dsigma2
#
# Version 0.1 (Apr 2021)
# ----------------------
#   Original code


dlogPhi.dmu <- function(y, mu, sigma, lower.tail=T)
{
  mu.sign <- ifelse(lower.tail, -1, 1)
  z <- (mu - y) / sigma * mu.sign
  # sum(varphi(z)$f) / sigma
  
  mu.sign * sum(varphi(z)$f) / sigma
} # end of dlogPhi.dmu


d2logPhi.dmu2 <- function(y, mu, sigma, lower.tail=T)
{
  mu.sign <- ifelse(lower.tail, -1, 1)
  
  z <- (mu - y)/sigma * mu.sign
  vphi <- varphi(z)
  
  sum(vphi$fp) / sigma^2
} # end of d2logPhi.dmu2


d2logPhi.dmudsigma <- function(y, mu, sigma, lower.tail=T)
{
  mu.sign <- ifelse(lower.tail, -1, 1)  

  z <- mu.sign * (mu - y) / sigma 
  vphi <- varphi(z)

  #mu.sign * sum(vphi$f + z*vphi$fp) / sigma^2
  
  - mu.sign * sum(z*vphi$fp + vphi$f) / sigma**2
} # end of d2logPhi.dmudsigma


dlogPhi.dsigma <- function(y, mu, sigma, lower.tail=T)
{
  mu.sign <- ifelse(lower.tail, -1, 1)
  z <- mu.sign * (mu - y) / sigma

  - sum(z * varphi(z)$f) / sigma
} # end of dlogPhi.dsigma


d2logPhi.dsigma2 <- function(y, mu, sigma, lower.tail=T)
{
  mu.sign <- ifelse(lower.tail, -1, 1)
  z <- mu.sign * (mu - y) / sigma
  vphi <- varphi(z)

  sum(vphi$fp*z^2 + 2*z*vphi$f) / sigma^2
} # end of d2logPhi.dsigma2


dlogPhiInterval.dmu <- function(gt, lt, mu, sigma)
{
  # gt, lt:       vectors of same length
  # mu, sigma:  scalars 

  z2 <- (lt-mu) / sigma
  z1 <- (gt-mu) / sigma

  num2.log <- dnorm(z2, log=T)
  num1.log <- dnorm(z1, log=T)
  num.log <- log.diff.exp(num1.log, num2.log)

  denom.log <- PhiInterval(gt, lt, mu, sigma, T)

  ratios <- exp(num.log$x - denom.log) * num.log$sign

  - sum(ratios) / sigma
} # end of dlogPhiInterval.dmu


d2logPhiInterval.dmu2 <- function(gt, lt, mu, sigma)
{
  z2 <- (lt-mu) / sigma
  z1 <- (gt-mu) / sigma
  log.sigma <- log(sigma)

  log.u <- log.dphi(z1, z2)
  log.v <- PhiInterval(z1, z2, 0, 1, T)

  log.u1 <- dnorm(z1, log=T)
  log.u2 <- dnorm(z2, log=T) 

  log.up1 <- list(x=log.u1 + log(abs(z1)) - log.sigma, sign=sign(z1))
  log.up2 <- list(x=log.u2 + log(abs(z2)) - log.sigma, sign=sign(z2))
    
  log.vp1 <- dPhi.dmu(gt, mu, sigma, log=T)
  log.vp2 <- dPhi.dmu(lt, mu, sigma, log=T)

  log.up <- log.diff(log.up1, log.up2)
  log.vp <- log.diff(log.vp1, log.vp2)
      
  log.d1 <- list(x=log.up$x + log.v, sign=log.up$sign)
  log.d2 <- list(x=log.u$x + log.vp$x, sign=log.u$sign * log.vp$sign)
  log.d <- log.diff(log.d2, log.d1)
  log.d$x <- log.d$x - 2 * log.v
  d <- exp(log.d$x) * log.d$sign

  - sum(d) / sigma
} # end of d2logPhiInterval.dmu2


d2logPhiInterval.dmudsigma <- function(gt, lt, mu, sigma)
{
  # let B = PhiInterval = Phi((b-mu)/sigma) - Phi((a-mu)/sigma)
  
  B <- PhiInterval(gt, lt, mu, sigma, T)

  m2 <- d2Phi.dmudsigma(lt, mu, sigma, log=T)
  m1 <- d2Phi.dmudsigma(gt, mu, sigma, log=T)
  d2B.dmudsigma <- log.diff(m1, m2) 
  
  u2 <- dPhi.dsigma(lt, mu, sigma, log=T)
  u1 <- dPhi.dsigma(gt, mu, sigma, log=T)
  Bp.sigma <- log.diff(u1, u2)

  v2 <- dPhi.dmu(lt, mu, sigma, log=T)
  v1 <- dPhi.dmu(gt, mu, sigma, log=T)
  Bp.mu <- log.diff(v1, v2)

  term2 <- log.mult(Bp.mu, Bp.sigma)
    term2$x <- term2$x - 2*B
    
  d2logPhiInterval.dmudsigma <- d2B.dmudsigma
    d2logPhiInterval.dmudsigma$x <- d2logPhiInterval.dmudsigma$x - B
    
   
  d2logPhiInterval.dmudsigma <- log.diff(term2, d2logPhiInterval.dmudsigma)

  sum(exp(d2logPhiInterval.dmudsigma$x) * d2logPhiInterval.dmudsigma$sign)
} # end of d2logPhiInterval.dmudsigma


dlogPhiInterval.dsigma <- function(gt, lt, mu, sigma)
{
  # gt, lt:       vectors of same length
  # mu, sigma:  scalars 

  u2 <- dPhi.dsigma(lt, mu, sigma, log=T)
  u1 <- dPhi.dsigma(gt, mu, sigma, log=T)
  u <- log.diff(u1, u2)

  log.v <- PhiInterval(gt, lt, mu, sigma, T)

  log.ratios <- u$x - log.v
  ratios <- u$sign * exp(log.ratios)

  sum(ratios)
} # end of dlogPhiInterval.dsigma


d2logPhiInterval.dsigma2 <- function(gt, lt, mu, sigma)
{
  z1 <- (gt - mu) / sigma
  z2 <- (lt - mu) / sigma
  
  B <- pnorm(z2) - pnorm(z1)
  
  Bp1 <- -z1 * dnorm(z1) / sigma
  Bp2 <- -z2 * dnorm(z2) / sigma
  Bp <- Bp2 - Bp1
  
  Bs1 <- -z1 * dnorm(z1) * (z1**2 - 2) / sigma**2
  Bs2 <- -z2 * dnorm(z2) * (z2**2 - 2) / sigma**2
  Bs <- Bs2 - Bs1
  
  d2sumLogB.dtheta2(B, Bp, Bs) 
} # end of d2logPhiInterval.dsigma2


d2phi.dmudsigma <- function(y, mu, sigma, log=F)
{
  z <- (y - mu) / sigma

  if (log)
  {
    b <- z^2 - 1
    x <- dnorm(z, log=T) + log(abs(b)) - 2 * log(sigma)
    s <- sign(b)
    out <- list(x=x, sign=s)
  }
  else
  {
    phi <- dnorm(z)
    out <- phi * (z^2 - 1) / sigma^2
  }

  out
} # end of d2phi.dmudsigma


dphi.dsigma <- function(y, mu, sigma, log=F)
{
  z <- (y - mu) / sigma

  if (log)
  {
    x <- dnorm(z, log=T) + 2*log(abs(z)) - log(sigma)
    out <- list(x=x, sign=rep(1,length(z)))
  }
  else
  {
    out <- dnorm(z) * z**2 / sigma
  }

  # if log = TRUE, return a list with dimensions (x, sign)
  # if log = FALSE, return a vector (or scalar, depending on whether y is a vector or a scalar)

  out
} # end of dphi.dsigma


d2phi.dsigma2 <- function(y, mu, sigma, log=F)
{
  z <- (y - mu) / sigma

  if (log)
  {
    log.phi <- dnorm(z, log=T)
    x <- log.phi + log(abs(z)) + log(abs(z^2-2)) - 2*log(sigma)
    out <- list(x=x, sign=sign(z) * sign(z^2 - 2))
  }
  else
  {
    phi <- dnorm(z)
    out <- phi * z *  (z^2 - 2) / sigma^2
  }

  out
} # end of d2phi.dsigma2

  
dPhi.dmu <- function(y, mu, sigma, lower.tail=T, log=F)
{
  mu.sign <- ifelse(lower.tail, -1, 1)
  
  z <- (mu - y)/sigma  * mu.sign
    

  if (log)
  {
    x <- dnorm(z, log=T) - log(sigma)
    out <- list(x=x, sign=rep(mu.sign, length(x)))
  }
  else
  {
    phi <- dnorm(z)
    out <- mu.sign * phi / sigma
  }
  
  out
} # end of dPhi.dmu


d2Phi.dmudsigma <- function(y, mu, sigma, lower.tail=T, log=F)
{
  mu.sign <- ifelse(lower.tail, -1, 1)
  z <- (mu - y) / sigma * mu.sign
  
  if (log)
  {
    f <- z**2 - 1
    l <- list(x=log(abs(f)), sign=sign(f))
    
    l2 <- dnorm(z, 0, 1, T) - 2*log(sigma)
    l$x <- l$x + l2
    l$sign <- l$sign * mu.sign
    l
  }
  else
  {
    mu.sign * dnorm(z) * (z**2-1) / sigma**2
  }
} # end of d2Phi.dmudsigma


dPhi.dsigma <- function(y, mu, sigma, lower.tail=T, log=F)
{
  mu.sign <- ifelse(lower.tail, -1, 1)
  
  a <- (mu - y) * mu.sign
  
  z <- a/sigma

  if (log)
  {
    x <- dnorm(z, log=T) + log(abs(a)) - 2*log(sigma)
    out <- list(x=x, sign=-sign(a))
  }
  else
  {
    phi <- dnorm(z)
    out <- - phi * z / sigma
  }

  out
} # end of dPhi.dsigma


d2Phi.dsigma2 <- function(y, mu, sigma, lower.tail=T, log=F)
{
  mu.sign <- ifelse(lower.tail, -1, 1)
  z <- (mu - y) / sigma * mu.sign
  
  if (log)
  {
   out <- dnorm(z, log=T) + log(abs(z)) + log(abs(z**2-2)) - 2*log(sigma)
   out <- list(x=out, sign=-sign(z)*sign(z**2-2))
  }
  else
  {
    phi <- dnorm(z)
    out <- - z * phi * (z**2 - 2) / sigma**2
  }
  
  out
} # end of d2Phi.dsigma2


dPhiInterval.dmu <- function(a, b, mu, sigma, return.as.log=F)
{
  if (return.as.log)
  {
    out <- phi.diff(a, b, mu, sigma, T)
    list(x=out$x-log(sigma), sign=-out$sign) 
  }
  else -1/sigma * phi.diff(a, b, mu, sigma)
} # end of dPhiInterval.dmu


d2PhiInterval.dmu2 <- function(a, b, mu, sigma, return.as.log=F)
{
  z2 <- (b - mu) / sigma
  z1 <- (a - mu) / sigma
  
  if (return.as.log)
  {
    phi1 <- -z1**2 / 2
    phi2 <- -z2**2 / 2 
  
    f1 <- list(x=phi1 + log(abs(z1)), sign=sign(z1))
    f2 <- list(x=phi2 + log(abs(z2)), sign=sign(z2))
    diff <- log.diff(f1, f2)
    list(x=diff$x - log(2*pi)/2 - 2*log(sigma), sign=-diff$sign)
  }
  else
  {
    phi2 = exp(-z2**2/2)
    phi1 = exp(-z1**2/2)
  
    - (z2*phi2 - z1*phi1) / sqrt(2*pi) / sigma**2
  }
} # end of d2PhiInterval.dmu2


dPhiInterval.dsigma <- function(a, b, mu, sigma, return.as.log=F)
{
  a <- dPhi.dsigma(a, mu, sigma, log=T)
  b <- dPhi.dsigma(b, mu, sigma, log=T)
  
  diff <- log.diff(a, b)
  
  if (return.as.log) diff
  else exp(diff$x) * diff$sign
} # end of dPhiInterval.dsigma


d2PhiInterval.dsigma2 <- function(a, b, mu, sigma, return.as.log=F)
{
  a <- d2Phi.dsigma2(a, mu, sigma, T, T)
  b <- d2Phi.dsigma2(b, mu, sigma, T, T) 
  
  diff <- log.diff(a, b)
  
  if (return.as.log) diff
  else exp(diff$x) * diff$sign
} # end of d2PhiInterval.dsigma2


dsumLogB.dtheta <- function(B, Bp, log=F)
{
  # If log=T, the two objects B & Bp consist in log-notation objects

  if (log) sum(exp(Bp$x - B$x)*Bp$sign*B$sign)
  else sum(Bp/B)
} # end of dsumLogB.dtheta


d2sumLogB.dtheta2 <- function(B, Bp, Bs, log=F)
{
  # If log=T, the objects B, Bp & Bs consist in log-notation objects
  
  if (log) sum(exp(Bs$x - B$x)*Bs$sign*B$sign) - sum(exp(2*(Bp$x-B$x)))
  else sum(Bs/B) - sum((Bp/B)**2)
} # end of d2sumLogB.dtheta2


d2sumLogB.dtheta1dtheta2 <- function(B, Bp1, Bp2, Bm)
{
  sum(Bm/B) - sum(Bp1*Bp2/B**2)
} # end of d2sumLogB.dtheta1dtheta2


# ----- Non-derivative useful functions ----------------------------------------


varphi <- function(z)
{
  log.phi <- dnorm(z, log=T)
  log.Phi <- pnorm(z, log.p=T)

  vphi <- exp(log.phi - log.Phi)
  vphi.prime <- - (z*vphi + vphi^2)

  list(f=vphi, fp=vphi.prime)
} # end of varphi