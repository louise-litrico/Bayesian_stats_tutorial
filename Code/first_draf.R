# Final assessment (tutorial) # 

# Libraries ----
library(tidyverse)
library(brms)
library(tidybayes)

library(ggeffects)
library(bayesplot)
library(LearnBayes)

# Load and check data ----
LPI_data <- read_csv("LPI_data.csv")
# head(LPI_data)  # 69 variables 
# str(LPI_data)  # in wide format now so will need to change that into tidy data

# Select subset and urn data into long format ----
France <- LPI_data %>% 
  gather(data = ., key = "year", value = "pop", 25:69) %>%  # turn into long format with one observation in each line
  filter(is.finite(pop)) %>%  # take out empty observations 
  mutate(year = parse_number(year)) %>%  # we need year as a numerical variable for our analysis
  filter(Country.list == "France", Common.Name == "Knot / Red knot") %>% 
  mutate(Location.of.population = case_when(grepl("Atlantic", Location.of.population) ~ "Atlantic Coast", 
                              grepl("Channel", Location.of.population) ~ "Channel Coast")) %>% 
  group_by(year) %>% 
  arrange(desc(year)) %>% 
  ungroup()

write.csv(France, "Data/red_knot.csv")

# length(unique(LPI_long$Country.list))  # lots of different locations included in this dataset, will need to refine 
# unique(LPI_long$Class)  # also lots of different species observed so will subset for that as well 
# make that into the original question we are trying to answer with the model not the other way around 

# unique(France$year)
# unique(France$Location.of.population)
# data from only 2 location but very big so might not need to account for that as a random effect
# = the whole of the west coast and the whole of the north coast of the country 
# maybe relevant since two different water bodies 

# Visualize the data distribution ---
# our population follows a poisson distribution because it is count data 
(hist_france <- ggplot(France, aes(x = pop)) +
    geom_histogram(colour = "#8B5A00", fill = "#CD8500") +
    theme_bw() +
    ylab("Count\n") + 
    xlab("\nCalidris canutus abundance") + 
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain")))              

# ggsave(filename = "histogram_redknot.png", hist_france, device = "png")

(boxplot_location <- ggplot(France, aes(Location.of.population, pop)) +
  geom_boxplot() +  # could be a significant effect between locations so should look at that 
  theme_bw() +
  xlab("Location\n") + 
  ylab("\nCalidris canutus abundance") + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "plain")))  

# ggsave(filename = "boxplot_location.png", boxplot_location, device = "png")

# First basic model ----
# How has knot populations changed over time in France? 
france1_mbrms <- brms::brm(pop ~ I(year - 1976), 
                         data = France, family = poisson(), chains = 3, 
                         iter = 3000, warmup = 1000)

summary(france1_mbrms)
fixef(france1_mbrms)  # to get more detailed values 

# use coef(model) if you have group-level effects (hierarchical data)
# explain what each term means (with graphs!)
# posterior mean for each effect 
# + lower and upper 95% credible intervals of the distribution for each effect 
# + effective sample size for each random effect 
# effective sample size should be high (1000 to 2000)
# credible intervals should not include 0 so we can say fixed effect is significant
# the narrower the interval, the more precise the estimate of the effect is 
# for random effects, we estimate the variance = effect is significant when variance histogram is not pushed against 0 
# the larger the spread, the less precise the estimation is 

plot(france1_mbrms)
pp_check(france1_mbrms)  # posterior predictive checks 

# Increasing the complexity of our model ----
# Are population counts related from year to year?
france2_mbrms <- brms::brm(pop ~ I(year - 1976) + (I(year - 1976)), 
                         data = France, family = poisson(), chains = 3, 
                         iter = 3000, warmup = 1000)

summary(france2_mbrms)
plot(france2_mbrms)

# Potential problems ----
# Trying to find out if population on the Atlantic coast are significantly more 
# abundant/increasing than those on the Channel coast?
# BUT the model doesn't work well 
# Solution 1 = change location into a numerical variable
France <- France %>% 
  mutate(location = case_when(grepl("Atlantic", Location.of.population) ~ 1, 
                              grepl("Channel", Location.of.population) ~ 2))

france3_mbrms <- brms::brm(pop ~ I(year - 1976) + location, 
                         data = France, family = poisson(), chains = 3, 
                         iter = 3000, warmup = 1000)

# Solution 2 = scale year to make the model work better (also doesn't work)
France$year.scaled <- scale(I(France$year - 1976), center = T)  # scaling time
France$pop.scaled <- scale(France$pop, center = T)  # scaling abundance

hist(France$pop.scaled)  # as you can see that the distribution changed 
# so will have to change it in the model as well
france4_mbrms <- brms::brm(pop.scaled ~ year.scaled + (1|location), 
                           data = France, family = gaussian(), chains = 3, 
                           iter = 3000, warmup = 1000)

# Solution 3 = setting a better prior
prior1 <- c(set_prior(prior = 'normal(0,6)', class='b', coef='year'), 	
            # global slope belongs to a normal distribution centered around 0
            set_prior(prior = 'normal(0,6)', class='b', coef='location'),
            # global slope for the second fixed effect
            set_prior(prior = 'normal(0,6)', class='Intercept', coef=''))  
            # global intercept
            # set_prior(prior = 'cauchy(0,2)', class='sd'))		# group-level intercepts and slopes

france5_mbrms <- brms::brm(pop ~ year + location, data = France, 
                           family = poisson(), chains = 3, prior = prior1,
                           iter = 3000, warmup = 1000)

# explain what the basic one looks like and what it does
# explain what a better one would be and how to define it 

# Solution 4 = increase iterations and warmup 

# Plotting the model results ----
loo(france1_mbrms,france3_mbrms, compare = TRUE)  # to assess which model fits the data better 
# look at elpd estimate for each model, the higher value the better the fit 
# so here our third model is better, so we sill carry on using that one 

summary(france3_mbrms)
plot(france3_mbrms)  # fuzzy caterpillar is present 
pp_check(france3_mbrms) # model fits the data very nicely! 

# model + raw data + CI
(model_fit <- France %>%
    add_predicted_draws(france3_mbrms) %>%
    ggplot(aes(x = year, y = pop)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),
                    alpha = 0.5, colour = "black") +
    geom_point(data = France, colour = "darkseagreen4", size = 3) +
    scale_fill_brewer(palette = "Greys") + 
    ylab("Calidris canutus abundance\n") +
    xlab("\nYear") +
    theme_bw() + 
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))

# ggsave(filename = "france3_fit.png", model_fit, device = "png")

# same as before but with a slope for each location
(location_fit <- France %>%
  group_by(Location.of.population) %>%
  add_predicted_draws(france3_mbrms) %>%
  ggplot(aes(x = year, y = pop, color = ordered(Location.of.population), fill = ordered(Location.of.population))) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1/4) +
  geom_point(data = France) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  ylab("Calidris canutus abundance\n") +
  xlab("\nYear") +
  theme_bw() + 
  theme(legend.title = element_blank()))

# ggsave(filename = "france3_location_fit.png", location_fit, device = "png")

# Trying out other brms stuff ----
# creating an estimate for a specific value but need to check if we can extrapolate the model
newdata <- data.frame(year = 2011, location = 1)
predict(france3_mbrms, newdata = newdata, re_formula = NA)
# same way to do the above but more precise because based on predictions of the regression line 
# rather than using a prediction function for each response
fitted(france3_mbrms, newdata = newdata, re_formula = NA)

# Setting a prior based on quantile information 
curve(dbeta(x,52.22,9.52105105105105)) # plot the prior
# The Beta distribution with parameters shape1 = a and shape2 = b has density
# Γ(a+b)/(Γ(a)Γ(b))x^(a-1)(1-x)^(b-1)
# for a > 0, b > 0 and 0 ≤ x ≤ 1 where the boundary values at x=0 or x=1 are defined as by continuity (as limits). 
# The mean is a/(a+b) and the variance is ab/((a+b)^2 (a+b+1)). 

france3_posterior <- france3_mbrms %>% 
  as_tibble() %>% 
  rename(year = "IyearM1976")

ggplot(france3_posterior, aes(x=b_IyearM1976)) + 
  geom_histogram()

# Watanabe-Akaike information criterion:
#   The WAIC has the advantages of:
#   
#   Averaging the likelihood over the posterior distribution rather than using the mean
# 
# Does not assume a multivariate Gaussian posterior distribution, as does the DIC (and AIC)
waic()  # we use loo but you can use this one as well 

success <- 0:25

plot(success, dpois(success, lambda=10),
     type='h',
     ylab='Probability',
     xlab ='Abundance',
     lwd=3)
