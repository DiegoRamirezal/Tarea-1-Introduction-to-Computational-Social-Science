library(tidyverse)
library(haven)
library(MatchIt)
library(Zelig)
library(marginaleffects)
library(cem)
library(estimatr)


##### calculo ATE directo####


read_data <- function(df)
{
  full_path <- paste("https://raw.github.com/scunning1975/mixtape/master/", 
                     df, sep = "")
  df <- read_dta(full_path)
  return(df)
}

nsw_dw <- read_data("nsw_mixtape.dta")

nsw_dw %>% 
  filter(treat == 1) %>% 
  summary(re78)

mean1 <- nsw_dw %>% 
  filter(treat == 1) %>% 
  pull(re78) %>% 
  mean()

nsw_dw$y1 <- mean1

nsw_dw %>% 
  filter(treat == 0) %>% 
  summary(re78)

mean0 <- nsw_dw %>% 
  filter(treat == 0) %>% 
  pull(re78) %>% 
  mean()

nsw_dw$y0 <- mean0

ate <- unique(nsw_dw$y1 - nsw_dw$y0)

nsw_dw <- nsw_dw %>% 
  filter(treat == 1) %>% 
  select(-y1, -y0)



#####Obtención puntaje de propensión e histogramas#######


nsw_dw_cpscontrol <- read_data("cps_mixtape.dta") %>% 
  bind_rows(nsw_dw) %>% 
  mutate(agesq = age^2,
         agecube = age^3,
         educsq = educ*educ,
         u74 = case_when(re74 == 0 ~ 1, TRUE ~ 0),
         u75 = case_when(re75 == 0 ~ 1, TRUE ~ 0),
         interaction1 = educ*re74,
         re74sq = re74^2,
         re75sq = re75^2,
         interaction2 = u74*hisp)

# estimating
logit_nsw <- glm(treat ~ age + agesq + agecube + educ + educsq + 
                   marr + nodegree + black + hisp + re74 + re75 + u74 +
                   u75 + interaction1, family = binomial(link = "logit"), 
                 data = nsw_dw_cpscontrol)

nsw_dw_cpscontrol <- nsw_dw_cpscontrol %>% 
  mutate(pscore = logit_nsw$fitted.values)

# mean pscore 
pscore_control <- nsw_dw_cpscontrol %>% 
  filter(treat == 0) %>% 
  pull(pscore) %>% 
  mean()

pscore_treated <- nsw_dw_cpscontrol %>% 
  filter(treat == 1) %>% 
  pull(pscore) %>% 
  mean()

# histogram
nsw_dw_cpscontrol %>% 
  filter(treat == 0) %>% 
  ggplot() +
  geom_histogram(aes(x = pscore), fill="blue" )

nsw_dw_cpscontrol %>% 
  filter(treat == 1) %>% 
  ggplot() +
  geom_histogram(aes(x = pscore), fill="blue",color="white")



  

##### Nearest-neighbor matching#####



m_out <- matchit(treat ~ age + agesq + agecube + educ +
                   educsq + marr + nodegree +
                   black + hisp + re74 + re75 + u74 + u75 + interaction1,
                 data = nsw_dw_cpscontrol, method = "nearest", 
                 distance = "logit", ratio =5)

m_data <- match.data(m_out)

z_out <- zelig(re78 ~ treat + age + agesq + agecube + educ +
                 educsq + marr + nodegree +
                 black + hisp + re74 + re75 + u74 + u75 + interaction1, 
               model = "ls", data = m_data)

x_out <- setx(z_out, treat = 0)
x1_out <- setx(z_out, treat = 1)

s_out <- sim(z_out, x = x_out, x1 = x1_out)

summary(s_out)



###### Nearest-neighbor matching, en base a codding de ayudantia ####

m_out <- matchit(treat ~ age + agesq + agecube + educ + educsq + marr + nodegree +
                   black + hisp + re74 + re75 + u74 + u75 + interaction1,
                 data = nsw_dw_cpscontrol, method = "nearest", 
                 distance = "logit", ratio = 5)

summary(m_out)

m_data <- match.data(m_out)

summary(m_data)


models <- list( regression_match_nocontrols <- lm(re78 ~ treat, data = m_data),
                regression_match_controls  <- lm(re78 ~ treat + age + agesq + 
                                                  agecube + educ + educsq + marr + nodegree +
                                                  black + hisp + re74 + re75 + u74 + u75 +interaction1,
                                                       weights=weights, data = m_data)
)

summary(regression_match_controls)
summary(regression_match_nocontrols)

comp  <- comparisons(regression_match_controls,
                     variables = "treat",
                     vcov = ~ subclass,
                     newdata = subset(m_data, treat == 1))



summary(comp)






#######Weighting on the propensity score#####




#continuation
N <- nrow(nsw_dw_cpscontrol)
#- Manual with non-normalized weights using all data
nsw_dw_cpscontrol <- nsw_dw_cpscontrol %>% 
  mutate(d1 = treat/pscore,
         d0 = (1-treat)/(1-pscore))

s1 <- sum(nsw_dw_cpscontrol$d1)
s0 <- sum(nsw_dw_cpscontrol$d0)


nsw_dw_cpscontrol <- nsw_dw_cpscontrol %>% 
  mutate(y1 = treat * re78/pscore,
         y0 = (1-treat) * re78/(1-pscore),
         ht = y1 - y0)

#- Manual with normalized weights
nsw_dw_cpscontrol <- nsw_dw_cpscontrol %>% 
  mutate(y1 = (treat*re78/pscore)/(s1/N),
         y0 = ((1-treat)*re78/(1-pscore))/(s0/N),
         norm = y1 - y0)

nsw_dw_cpscontrol %>% 
  pull(ht) %>% 
  mean()

nsw_dw_cpscontrol %>% 
  pull(norm) %>% 
  mean()

#-- trimming propensity score
nsw_dw_cpscontrol <- nsw_dw_cpscontrol %>% 
  select(-d1, -d0, -y1, -y0, -ht, -norm) %>% 
  filter(!(pscore >= 0.9)) %>% 
  filter(!(pscore <= 0.1))

N <- nrow(nsw_dw_cpscontrol)

#- Manual with non-normalized weights using trimmed data
nsw_dw_cpscontrol <- nsw_dw_cpscontrol %>% 
  mutate(d1 = treat/pscore,
         d0 = (1-treat)/(1-pscore))

s1 <- sum(nsw_dw_cpscontrol$d1)
s0 <- sum(nsw_dw_cpscontrol$d0)

nsw_dw_cpscontrol <- nsw_dw_cpscontrol %>% 
  mutate(y1 = treat * re78/pscore,
         y0 = (1-treat) * re78/(1-pscore),
         ht = y1 - y0)

#- Manual with normalized weights with trimmed data
nsw_dw_cpscontrol <- nsw_dw_cpscontrol %>% 
  mutate(y1 = (treat*re78/pscore)/(s1/N),
         y0 = ((1-treat)*re78/(1-pscore))/(s0/N),
         norm = y1 - y0)

nsw_dw_cpscontrol %>% 
  pull(ht) %>% 
  mean()

nsw_dw_cpscontrol %>% 
  pull(norm) %>% 
  mean()





######## Coarsened exact matching#####



m_out <- matchit(treat ~ age + agesq + agecube + educ +
                   educsq + marr + nodegree +
                   black + hisp + re74 + re75 + 
                   u74 + u75 + interaction1,
                 data = nsw_dw_cpscontrol, 
                 method = "cem", 
                 distance = "logit")

m_data <- match.data(m_out)

m_ate <- lm_robust(re78 ~ treat, 
                   data = m_data,
                   weights = m_data$weights)


summary(m_ate)

