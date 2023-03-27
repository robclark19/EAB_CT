#### SET UP ####

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(ggeffects)
library(showtext)

font_add_google(name="Open Sans", family="Open Sans")
showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

eab = read_xlsx("./Data/Extra/CleanPropdata.xlsx",sheet="EABCountR") |>
  gather(key="year",value="count_eab",`0`:`10`)
  
prey = read_xlsx("./Data/Extra/CleanPropdata.xlsx",sheet="AllCountR") |>
  gather(key="year",value="count_all",`0`:`10`)

dat = left_join(eab,prey,by=c('Site','year')) |>
  mutate(
    prop_eab = count_eab/count_all,
    year = as.double(year)
  ) |>
  as.data.frame()


#### EAB INCIDENCE VS TIME ####
# Full binomial glmer: quadratic term, random slope and intercepts
full_model = glmer(prop_eab ~ year + I(year^2) + (1+year+I(year^2)|Site),
                  data=dat, weights=count_all,
                  family=binomial)
summary(full_model)

# Compare against simpler models
m2 = update(full_model, .~ year + I(year^2) + (1+year|Site))
m1 = update(full_model, .~ year + I(year^2) + (year|Site))
m0 = update(full_model, .~ year + I(year^2) + (1|Site))
anova(full_model,m2,m1,m0)
# full model favoured by a wide margin

# find the x value of the vertex (-b/2a)
a = summary(full_model)$coefficients[,1][3]
b = summary(full_model)$coefficients[,1][2]
-b/(2*a)
# peak is at year = 4.05

# Plot the parabola
basesize = 7
ggpredict(full_model,terms="year [all]") %>%
  ggplot(aes(x=x,y=predicted)) +
  geom_line(size=0.2) +
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.1) +
  labs(title="",y="Proportion of EAB in Prey",x="Year Post-Detection") +
  scale_x_continuous(breaks=seq(0,10,by=2)) +
  theme_bw() + theme(
    text = element_text(size=basesize, family="Open Sans"),
    axis.title = element_text(size=1.2*basesize),
    panel.grid = element_blank(),
    axis.ticks = element_line(size=0.1)
  )
ggsave("./R/outputs/figure2.png",units="px",dpi=300,width=800,height=500,device=ragg::agg_png())

# Plot at the site level
ggpredict(full_model, terms=c("year","Site[all]"), type="re") |>
  plot(ci=FALSE, limit.range=TRUE) + see::scale_color_flat()


#### Speed of invasion ####
library(invasionSpeed)

speed = read.csv("./Data/Extra/speed_data.csv") %>%
  filter(!is.na(FirstYearAny))

# Specify coordinates and years
coord = cbind(speed$Long,speed$Lat)
dates = speed$FirstYearAny

# Set default MCMC parameters
amcmc = list(n.batch=10,batch.length=100,accept.rate=0.3)

# Obtain localGrad
# For Albers specify standard parallels containing the study area (State of CT)
out.grad2 = localgrad(dates=dates, coord.longlat=coord, Albers=c(40.5,42.5), n.samp=1000, amcmc=amcmc)
summary.localgrad(out.grad2) #EAB

par(mar=c(5,6,4,1)+.1)
plotgrad(out.grad2,
         cex.axis=0.5,cex.lab=0.75,
         pch=".",database="state",
         main="",
         xlim=c(-74.5,-71),ylim=c(40.5,42.5))
#Mean speed of spread is 20.53 km/yr
#Median speed of spread if 17.01 km/yr

#Standard error of spread:
means.km = out.grad2$conv.factor*out.grad2$means
mag.speed = sqrt( means.km[,1]^2 + means.km[,2]^2 )
se = sd(mag.speed)/sqrt(length(mag.speed))
se #2.09 km/yr
