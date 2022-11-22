## 'Capstone project 2: Choose Your Own' - Alicia Alfonso ##

# I. Introduction

# Daphniids are small crustaceans that live in freshwater ecosystems. Like
# many aquatic invertebrates, their body is composed of calcium (Ca) and
# their survival likely depends on the availability of this element in the
# water. The dataset used for this project are simplified results from a
# personal experiment (Alfonso, A. 2018 unpublished). The aim of the
# experiment was to assess the inter-clonal variation between two
# populations/clones (Northern and Southern) at different Ca
# concentrations (Ca 0.5 mg/l, 4 mg/l and 106 mg/l) and temperatures (15
# and 25ºC). With this information, we will create a model to predict the
# survival of Daphniids evaluating the effect of the different variables
# measured.


# II. Methods and analysis

## Part 1: Data exploration
data <- read.table("CapstoneProjectData.txt", header=TRUE) # Load dataset
library(survival) # For modelling
library(survminer) # For survival curves
library(dplyr) # For analysis
library(ggplot2) # For visualization

str(data) # 123 obs. of  6 variables
head(data, 10)
summary(data)


p <- data %>%
  ggplot(aes(x=Time)) +
  geom_histogram( binwidth=3, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Final age") +
  theme(plot.title = element_text(size=15))
p

sum(data$Time == 88) # 13 alive at end

data %>%
  group_by(Clone, Temp, Ca) %>%
  summarise_at(vars(Time), list(name = mean))


## Part 2: Modelling survival
s <- Surv(data$Time, data$Dead) # Creating our hazard rate
# Survival Plots
m1 <- survfit(s ~ 1, data = data, conf.type = "log-log") 
m1p <- ggsurvplot(m1, # Our base curve for all treatments
                  data = data,
                  palette = "#2E9FDF",
                  title = "Total survival",
                  size = .5,
                  ggtheme = theme_bw())

m2 <- survfit(s ~ Clone, data = data, conf.type = "log-log")
m2p <- ggsurvplot(m2, # Our curve per clonal population
                  data = data,
                  conf.int=TRUE,
                  title = "Survival per Clone",
                  size = .5,  # change line size
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  palette = c("#A569BD", "#48C9B0"), # custom color palette
                  pval = TRUE)

m3 <- survfit(s ~ Temp, data = data, conf.type = "log-log")
m3p <- ggsurvplot(m3, # Our curve by temperature
                  data = data,
                  conf.int=TRUE,
                  title = "Survival per Temperature",
                  size = .5,  # change line size
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  palette = c("#006399", "#CC0000"), # custom color palette
                  pval = TRUE)

m4 <- survfit(s ~ Ca, data = data, conf.type = "log-log")
m4p <- ggsurvplot(m4, # Our curve per Ca concentration
                  data = data,
                  conf.int=TRUE,
                  title = "Survival per Ca concentration",
                  size = .5,  # change line size
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  palette = c("#F39C12","#E67E22","#D35400"),
                  pval = TRUE)

survplots <- list(m1p, m2p, m3p, m4p)
arrange_ggsurvplots(survplots, # Arrange multiple survplots
                    ncol = 2, nrow = 2)
# Median survival time is 39 days
# Overlapping between clones
# Overlapping between temperatures
# No-overlapping for 0.5 mg Ca/l


# IV. Results

data$Clone <- as.factor(data$Clone)
data$Temp <- as.factor(data$Temp)
data$Ca <- as.factor(data$Ca)

# Non-parametric Cox proportional hazards models
# Fits a Cox proportional hazards regression model
summary(mod1 <- coxph(s ~ Clone + Temp * Ca, data=data)) # Ca:Temp interaction
summary(mod2 <- coxph(s ~ Clone + Temp + Ca, data=data)) # No interaction
summary(mod3 <- coxph(s ~ (Clone + Temp + Ca)^2, data=data)) # Ca:Temp, Clone:Temp, Clone:Ca interaction

anova(mod1, mod2, mod3) # (!) Model 3 is already significant
AIC(mod1, mod2, mod3)   # 

# Using the step method to refine our best model (mod3) into a better one (mod4)
summary(mod4 <- step(mod3))
# Lowest AIC: Clone:Temp & Clone:Ca retained, no Temp:Ca interaction
AIC(mod1, mod2, mod3, mod4)
# Creating the better fitted model adding interactions Clone x Temp and Clone x Ca
mod4 <- coxph(s ~ Clone + Temp + Ca + Clone:Temp + Clone:Ca, data=data)
summary(mod4)

# Model 4 visualization separated by treatment 
mod4S15 <- survfit(mod4, newdata=data.frame(Clone="Southern", Temp="15", Ca=c("0.5", "4", "106")))
pmod4S15 <- ggsurvplot(mod4S15, conf.int=TRUE,
                       title = "Survival Southern 15C",
                       legend.labs = c("0.5", "4", "106"),
                       data=data.frame(Clone="Southern", Temp="15", Ca=c("0.5", "4", "106")),
                       size = .5,  # change line size
                       xlab = "Time (days)",
                       ggtheme = theme_classic(), # Change ggplot2 theme
                       palette = c("#F39C12","#E67E22","#D35400"))

mod4S25 <- survfit(mod4, newdata=data.frame(Clone="Southern",Temp="25", Ca=c("0.5", "4", "106")))
pmod4S25 <- ggsurvplot(mod4S25, conf.int=TRUE,
                       title = "Survival Southern 25C",
                       legend.labs = c("0.5", "4", "106"),
                       data=data.frame(Clone="Southern", Temp="25", Ca=c("0.5", "4", "106")),
                       size = .5,  # change line size
                       xlab = "Time (days)",
                       ggtheme = theme_classic(), # Change ggplot2 theme
                       palette = c("#F39C12","#E67E22","#D35400"))

mod4N15 <- survfit(mod4, newdata=data.frame(Clone="Northern",Temp="15", Ca=c("0.5", "4", "106")))
pmod4N15 <- ggsurvplot(mod4N15, conf.int=TRUE,
                       title = "Survival Northern 15C",
                       legend.labs = c("0.5", "4", "106"),
                       data=data.frame(Clone="Northern", Temp="15", Ca=c("0.5", "4", "106")),
                       size = .5,  # change line size
                       xlab = "Time (days)",
                       ggtheme = theme_classic(), # Change ggplot2 theme
                       palette = c("#F39C12","#E67E22","#D35400"))

mod4N25 <- survfit(mod4, newdata=data.frame(Clone="Northern",Temp="25", Ca=c("0.5", "4", "106")))
pmod4N25 <- ggsurvplot(mod4N25, conf.int=TRUE,
                       title = "Survival Northern 25C",
                       legend.labs = c("0.5", "4", "106"),
                       data=data.frame(Clone="Northern", Temp="25", Ca=c("0.5", "4", "106")),
                       size = .5,  # change line size
                       xlab = "Time (days)",
                       ggtheme = theme_classic(), # Change ggplot2 theme
                       palette = c("#F39C12","#E67E22","#D35400"))

mod4plots <- list(pmod4S15, pmod4S25, pmod4N15, pmod4N25)
arrange_ggsurvplots(mod4plots, # Arrange multiple survplots
                    title = "Model 4",
                    ncol = 2, nrow = 2)

# V. Conclusion
# Modelling survival of organisms can give us a deep insight into which
# factors contribute to their success in an ecosystem. Using R for this
# purpose proves extremely useful and practical as we are able to generate
# predictions from smaller datasets. In our case, we establish the strong
# negative effect that low calcium concentrations have in the survival of
# Daphniids. We can conclude that the mortality will highly increase if
# there is a depletion of this element and that northern populations will
# be the most affected. Our model also hints towards strong interactions
# between variables that leads to more interesting questions regarding
# phenotypic plasticity or even toxicity for even higher calcium
# concentrations.
