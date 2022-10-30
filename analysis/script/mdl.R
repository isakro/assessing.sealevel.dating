library(tidyverse)
library(MASS)
library(poweRlaw)
library(AICcmodavg)

# Load simulation results
load(here::here("analysis/data/derived_data/06data.RData"))

# Vertical distances, including negative values
vertdistneg <- distances %>%
  filter(year <= -2500 & rcarb_cor == "t")
vertdistneg$vertdist <- vertdistneg$vertdist + 0.001

vertdisto <- distances %>%
  filter(year <= -2500 & rcarb_cor == "t") %>%
  mutate(vertdist = ifelse(vertdist < 0, 0, vertdist))

# Vertical distances, forcing negative values to zero
vertdist <- distances %>%
  filter(year <= -2500 & rcarb_cor == "t") %>%
  mutate(vertdist = ifelse(vertdist < 0, 0, vertdist))
vertdist$vertdist <- vertdist$vertdist + 0.001

# Vertical distances, excluding negative values
vertdistna <- distances %>%
  filter(year <= -2500 & rcarb_cor == "t") %>%
  mutate(vertdist = ifelse(vertdist < 0, NA, vertdist)) %>%
  filter(!is.na(vertdist))
vertdistna$vertdist <- vertdistna$vertdist + 0.001

# Fit normal distribution to all values
normfit <- MASS::fitdistr(vertdistneg$vertdist, "normal")
normdat <- data.frame(
  x = vertdistneg$vertdist,
  px = dnorm(vertdistneg$vertdist,
             mean =  normfit$estimate[1],
             sd =  normfit$estimate[1]))

# Fit exponential function to only positive values
expfit <- fitdistr(vertdist$vertdist, "exponential")
expdat <- data.frame(
  x = vertdist$vertdist,
  px = dexp(as.numeric(vertdist$vertdist), rate = expfit$estimate))
expdatcum <- data.frame(
  x = vertdist$vertdist,
  px = pexp(as.numeric(vertdist$vertdist), rate = expfit$estimate))

# Fit exponential function to only positive values
expfitna <- fitdistr(vertdistna$vertdist, "exponential")
expdatna <- data.frame(
  x = vertdistna$vertdist,
  px = dexp(as.numeric(vertdistna$vertdist), rate = expfitna$estimate))

# Fit log-normal function to only positive values
lognormfit <- fitdistr(vertdist$vertdist, "log-normal")
lognormdat <- data.frame(
  x = vertdist$vertdist,
  px = dlnorm(vertdist$vertdist, meanlog = lognormfit$estimate[1],
              sdlog = lognormfit$estimate[2]))

# Fit a continuous power law (using the poweRlaw package)
pl <- conpl$new(vertdist$vertdist)
pl$setXmin(min(vertdist$vertdist))
pl$setPars(estimate_pars(pl))
plot(pl)
pl_xy <- lines(pl)
pl_xy$px <- dist_pdf(pl, q = pl_xy$x)

# Fit gamma distrbution to only positive values. (Warnings are thrown as optim
# tries negative paramter values, but the parameters for gamma are
# positive real numbers).
gammafit <- fitdistr(vertdist$vertdist, "gamma")
gammadat <- data.frame(
  x = vertdist$vertdist,
  px = dgamma(vertdist$vertdist, shape = gammafit$estimate[1],
              rate = gammafit$estimate[2]))

# Fit poisson distrbution to only positive values
poissonfit <- fitdistr(vertdist$vertdist, "Poisson")
poissondat <- data.frame(
  x = vertdisto$vertdist,
  px = dgamma(vertdist$vertdist, shape = poissonfit$estimate[1]))

# Fit logistic distrbution to only positive values
logisticfit <- fitdistr(vertdist$vertdist, "logistic")
logisticdat <- data.frame(
  x = vertdist$vertdist,
  px = dlogis(vertdist$vertdist, location = logisticfit$estimate[1],
              scale = logisticfit$estimate[2]))

vertdistneg %>%
ggplot() + geom_histogram(aes(vertdist, y = ..density..), binwidth = 1,
                          fill = NA,
                          colour = "darkgrey") +
  geom_line(data = expdat, aes(x = x, y = px), col = "red") +
  geom_line(data = expdatna, aes(x = x, y = px), col = "gold1") +
  geom_line(data = normdat, aes(x = x, y = px), col = "black") +
  geom_line(data = pl_xy, aes(x = x, y = px), col = "darkgreen") +
  geom_line(data = lognormdat, aes(x = x, y = px), col = "darkred") +
  geom_line(data = gammadat, aes(x = x, y = px), col = "blue") +
  geom_line(data = logisticdat, aes(x = x, y = px), col = "purple") +
  ylim(c(0, 0.25)) +
  # geom_text(aes(x = Inf, y = Inf, hjust = 1.2,
  #               vjust = 5,
  #               label = label_text), col = "red", parse = TRUE) +
  theme_bw() +
  labs(title= "Vertical distance (m) \nBinwidth = 1m", x = "", y = "Density")

vertdist <- vertdist %>%
  mutate(rank = min_rank(vertdist),
         ccdf = round((length(rank)-rank+1)/length(rank), 3))

ggplot() +
geom_point(data = vertdist, aes(x = vertdist, y = cumsum(vertdist)))


vertdist %>%
  ggplot() + geom_density(aes(vertdist),
                            fill = NA,
                            colour = "darkgrey") +
  geom_line(data = expdat, aes(x = x, y = px), col = "red") +
  geom_line(data = expdatna, aes(x = x, y = px), col = "gold1") +
  geom_line(data = pl_xy, aes(x = x, y = px), col = "darkgreen") +
  geom_line(data = lognormdat, aes(x = x, y = px), col = "darkred") +
  geom_line(data = gammadat, aes(x = x, y = px), col = "blue") +
  geom_line(data = logisticdat, aes(x = x, y = px), col = "purple") +
  geom_line(data = poissondat, aes(x = x, y = px), col = "purple") +
  ylim(0, 5) +
  # geom_text(aes(x = Inf, y = Inf, hjust = 1.2,
  #               vjust = 5,
  #               label = label_text), col = "red", parse = TRUE) +
  theme_bw() +
  labs(title= "Vertical distance (m) \nBinwidth = 1m", x = "", y = "Density")


loglik <- c(logLik(normfit)[1], logLik(expfit)[1],
            logLik(lognormfit)[1], dist_ll(pl), logLik(gammafit)[1],
            logLik(logisticfit)[1])
dfs <- c(2, 1, 2, 1, 2, 2)
names <- c("Normal", "Exponential", "Log-normal",
           "Power law", "Gamma", "Logistic")

BIC <- bictabCustom(logL = loglik, K = dfs, modnames = names,
                    nobs = nrow(vertdist))

BIC

AIC <- aictabCustom(logL = log_lik, K = no_pars, modnames = modnames,
                    nobs = nobs)
AIC

