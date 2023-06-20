library(rstan)
library(tidyverse)
library(lubridate)

# Reading data and identifying the dates
df.pre  <- read_csv("data/building6pre.csv")  %>% mutate(Date = mdy_hm(Date))
df.dur  <- read_csv("data/building6during.csv")  %>% mutate(Date = mdy_hm(Date))
df.post <- read_csv("data/building6post.csv") %>% mutate(Date = mdy_hm(Date))

# Quick preprocessing and daily averaging
daily.average <- function(df) {
  df %>% 
    mutate(day = as_date(Date)) %>%
    group_by(day) %>% 
    summarise(T = (mean(OAT)-32)/1.8,
              E = sum(`Building 6 kW`),
              .groups = 'drop') %>% 
    mutate(wday = wday(day),
           week.end = wday==1 | wday==7)
}

# Filtering out the week ends
df.pre.day  <- daily.average(df.pre) %>% filter(!week.end) %>% mutate(period = 'pre')
df.dur.day  <- daily.average(df.dur) %>% filter(!week.end) %>% mutate(period = 'during')
df.post.day <- daily.average(df.post) %>% filter(!week.end) %>% mutate(period = 'post')

model_data <- list(
  N_pre = nrow(df.pre.day),
  x_pre = df.pre.day$T,
  y_pre = df.pre.day$E,
  N_post = nrow(df.post.day),
  x_post = df.post.day$T,
  y_post = df.post.day$E
)

# Fitting the model
fit <- stan(
  file = 'models/changepoint_ma15.stan',  # Stan program
  data = model_data,        # named list of data
  chains = 4,               # number of Markov chains
  warmup = 3000,            # number of warmup iterations per chain
  iter = 4000,              # total number of iterations per chain
  cores = 2,                # number of cores (could use one per chain)
  refresh = 0,              # progress not shown
  control = list(max_treedepth = 12),
  seed = 5678,
)

print(fit, pars = c("alpha", "beta_h", "tau_h", "beta_c", "tau_c", "sigma", "theta1", "theta5", "savings"))
#traceplot(fit, pars = c("alpha", "beta_h", "tau_h", "beta_c", "theta1", "theta5"))
#pairs(fit, pars = c("alpha", "beta_h", "tau_h", "beta_c", "tau_c", "sigma", "theta1", "theta5"))

la <- rstan::extract(fit)
predict.mean <- colMeans(la$y_pre_pred)

residuals <- predict.mean - df.pre.day$E      # residuals
ssres <- sum(residuals^2)                     # sum of squared residuals
sstot <- sum((df.pre.day$E-mean(df.pre.day$E))^2) # total sum of squares

R2 <- 1 - ssres / sstot   # coefficient of determination
CVRMSE <- sqrt(ssres / nrow(df.pre.day)) / mean(df.pre.day$E)

print(paste("R2 =", R2))
print(paste("CV(RMSE) =", CVRMSE))

quantile(la$savings, probs=c(0.025, 0.1, 0.5, 0.9, 0.975))

acf(residuals)

ggplot() +
  geom_point(mapping = aes(x=residuals, y=lag(residuals, n=1, default=0)))

ggplot() +
  geom_point(mapping = aes(x=residuals, y=lag(residuals, n=5, default=0)))