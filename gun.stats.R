# gun stats
# Corey Bradshaw
# October 2017
# data from https://theconversation.com/gun-control-in-america-by-the-right-and-wrong-numbers-49573

## Remove everything
rm(list = ls())

# Set functions
AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]}
  return(AICc.vec)}

k.glm <- function(x) {
  if (length(x$df.residual) == 0) k <- sum(x$dims$ncol) else k <- (length(x$coeff)+1)}

delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
chdev.glm <- function(x) ((( as.numeric(x[12]) - as.numeric(x[10]))/ as.numeric(x[12]))*100) ## % change in deviance, where x is glm object

setwd("~/Documents/Statistics/guns/") # set appropriately
gun.dat <- read.table("gun.csv",header=T,sep=",") # gun data

plot(gun.dat$gunpercap, gun.dat$deathinjp100k, pch=19)
plot(log10(gun.dat$gunpercap), log10(gun.dat$deathinjp100k), pch=19, xlab="log per capita gun ownership", ylab = "log death+injury/100K", xlim=c(0,2), ylim=c(0,2))
fit.pow0 <- lm(log10(gun.dat$deathinjp100k) ~ 0 + log10(gun.dat$gunpercap)) # force through origin y = bx
summary(fit.pow0)
fit.pow1 <- lm(log10(gun.dat$deathinjp100k) ~ log10(gun.dat$gunpercap)) # y = a + bx
summary(fit.pow1)
abline(fit.pow0)

AICc.weights <- weight.IC(delta.IC(AICc(fit.pow1,fit.pow0))) # AICc for each model
AICc.weights
ER <- AICc.weights[1]/AICc.weights[2] # information-theoretic evidence ratio of origin model (bx) vs. a+bx model
ER
