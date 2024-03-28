# boxplots for weight changes
# Jenn Kim

# load libraries
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(nlme)
library(lsmeans)
library(plotrix)

# load data
data <- read.csv("BehaviorCohort_BodyWeights.csv")
data <- as.data.frame(data)

data$LPS.Treatment <- factor(data$LPS.Treatment, levels=c("PBS","1XLPS","2XLPS","4XLPS"))
data$Sex <- factor(data$Sex) 

# calculate weight change relative to before experiment for each day
data2 <- data %>% 
  mutate(PercentWeightLoss_BeforeExperiment = round(((Before.Experiment - Before.Experiment)/Before.Experiment)*100,1)) %>%
  mutate(PercentWeightLoss_Day1 = round(((Before.Experiment - Day.1)/Before.Experiment)*100,1)) %>%
  mutate(PercentWeightLoss_Day2 = round(((Before.Experiment - Day.2)/Before.Experiment)*100,1)) %>%
  mutate(PercentWeightLoss_Day3 = round(((Before.Experiment - Day.3)/Before.Experiment)*100,1)) %>%
  mutate(PercentWeightLoss_Day4 = round(((Before.Experiment - Day.4)/Before.Experiment)*100,1)) %>%
  mutate(PercentWeightLoss_Day5 = round(((Before.Experiment - Day.5)/Before.Experiment)*100,1)) %>%
  mutate(PercentWeightLoss_Day6 = round(((Before.Experiment - Day.6)/Before.Experiment)*100,1)) %>%
  mutate(PercentWeightLoss_Day7 = round(((Before.Experiment - Day.7)/Before.Experiment)*100,1)) %>%
  mutate(PercentWeightLoss_Day8 = round(((Before.Experiment - Day.8)/Before.Experiment)*100,1)) %>%
  mutate(across(PercentWeightLoss_Day1:PercentWeightLoss_Day8, ~.*-1)) # for graphing purposes, * -1 makes weight loss negative and weight gain positive
  

# line graphs for body weight changes
data3 <- data2 %>% 
  gather(Timepoint, PercentWeightLoss, 12:ncol(data2)) %>% 
  group_by(LPS.Treatment, Timepoint) %>%
  summarise(mean = mean(PercentWeightLoss), sem = std.error(PercentWeightLoss))

# clean up day names
data3 <- data3 %>% mutate(Timepoint = 
                           case_when(Timepoint=="PercentWeightLoss_BeforeExperiment" ~ "Before Start",
                                     Timepoint=="PercentWeightLoss_Day1" ~ "Day 1",
                                     Timepoint=="PercentWeightLoss_Day2" ~ "Day 2",
                                     Timepoint=="PercentWeightLoss_Day3" ~ "Day 3",
                                     Timepoint=="PercentWeightLoss_Day4" ~ "Day 4",
                                     Timepoint=="PercentWeightLoss_Day5" ~ "Day 5",
                                     Timepoint=="PercentWeightLoss_Day6" ~ "Day 6",
                                     Timepoint=="PercentWeightLoss_Day7" ~ "Day 7",
                                     Timepoint=="PercentWeightLoss_Day8" ~ "Day 8",
                                     Timepoint=="PercentWeightLoss_Day9" ~ "Day 9"))

# run stats on individual values
library(lme4)
stats <- data2 %>% 
  rownames_to_column(var="AnimalID") %>%
  gather(Timepoint, PercentWeightLoss, 13:ncol(.))

stats$LPS.Treatment <- factor(stats$LPS.Treatment)
stats$Timepoint <- factor(stats$Timepoint)

# model
model <- lme(PercentWeightLoss ~ LPS.Treatment*Timepoint, ~1|AnimalID, data=stats)

# anova
anova = anova(model)

# posthoc 1: differences in timepoints within treatment
refgrid <- ref.grid(model)
ph <- emmeans(refgrid, ~Timepoint|LPS.Treatment)
ph2 <- contrast(ph, method="pairwise", adjust="none")
ph3 <- test(ph2, by=NULL, adjust="BH")
posthoc1 <- as.data.frame(ph3)

# posthoc 2: differences in treatment (only considering comparisons against PBS) within timepoint
stats <- data2 %>% 
  rownames_to_column(var="AnimalID") %>%
  gather(Timepoint, PercentWeightLoss, 13:ncol(.)) %>%
  mutate(LPS.Treatment = 
           case_when(LPS.Treatment=="1XLPS" ~ "LPSx1",
                     LPS.Treatment=="2XLPS" ~ "LPSx2",
                     LPS.Treatment=="3XLPS" ~ "LPSx3",
                     LPS.Treatment=="4XLPS" ~ "LPSx4",
                     LPS.Treatment=="PBS" ~ "PBS"))
stats$LPS.Treatment <- factor(stats$LPS.Treatment, levels=c("PBS","LPSx1","LPSx2","LPSx3","LPSx4"))
stats$Timepoint <- factor(stats$Timepoint)

model <- lme(PercentWeightLoss ~ LPS.Treatment*Timepoint, ~1|AnimalID, data=stats)

refgrid <- ref.grid(model)
ph <- emmeans(refgrid, ~LPS.Treatment|Timepoint)

PBS = c(1,0,0,0)
LPSx1 = c(0,1,0,0)
LPSx2 = c(0,0,1,0)
LPSx4 = c(0,0,0,1)

ph2 <- contrast(ph, 
                method=list("PBS - 1XLPS" = PBS - LPSx1,
                            "PBS - 2XLPS" = PBS - LPSx2,
                            "PBS - 4XLPS" = PBS - LPSx4), 
                adjust="none")
ph3 <- test(ph2, by=NULL, adjust="BH")
posthoc2 <- as.data.frame(ph3)

anova$Significant <- ifelse(anova$`p-value` < 0.05, "significant", "ns")
posthoc1$Significant <- ifelse(posthoc1$p.value < 0.05, "significant", "ns")
posthoc2$Significant <- ifelse(posthoc2$p.value < 0.05, "significant", "ns")

# save output
write.csv(anova, "WeightChanges_anova.csv")
write.csv(posthoc1, "WeightChanges_posthoc_WithinTreatment.csv")
write.csv(posthoc2, "WeightChanges_posthoc_WithinTimepointvsPBSonly.csv")

# plot with mean and sem as shading
pdf("BehaviorCohort_BodyWeightChange_Percentage2.pdf", width=7, height=5)
ggplot(data3, aes(x=Timepoint, y=mean, group=LPS.Treatment, color=LPS.Treatment)) +
  geom_line(linewidth=.5) +
  scale_color_viridis_d() +
  theme_cowplot(font_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem, fill=LPS.Treatment),alpha = .4) +
  scale_fill_viridis_d() +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "% weight change") +
  labs(fill="LPS Treatment", color="LPS Treatment") +
  ggtitle("Body weight changes 24 hrs after each experiment day")
dev.off()