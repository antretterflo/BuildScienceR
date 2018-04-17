# Adaptive Comfort standards ----------------------------------------------

# Functions to compute hygrothermal values
# Florian Antretter
# Date: 2018-04-12
# Version: 0.1

calcLimitTemp <- function(Tout, Slope = 0.33, Intercept = 18.8) {
  Slope * Tout + Intercept
}

calcLimitTemp(20)
library(ggplot2)

ggplot() +
  coord_cartesian(xlim = c(5, 30), ylim = c(15, 30)) +
  geom_path(aes(x = c(10, 30), y = calcLimitTemp(c(10, 30), 0.31, 21.3)), color = 'red') +
  geom_path(aes(x = c(10, 30), y = calcLimitTemp(c(10, 30), 0.31, 20.3)), color = 'red2', size = 1.3) +
  geom_path(aes(x = c(10, 30), y = calcLimitTemp(c(10, 30), 0.31, 17.8)), color = 'black', size = 1.3, linetype = 2) +
  geom_path(aes(x = c(10, 30), y = calcLimitTemp(c(10, 30), 0.31, 15.3)), color = 'blue2', size = 1.3) +
  geom_path(aes(x = c(10, 30), y = calcLimitTemp(c(10, 30), 0.31, 14.3)), color = 'blue') +
  labs(title = "ASHRAE 55", x = "Prevailing Outdoor Temperature [°C]", y = "Indoor Operative Temperature [°C]")

