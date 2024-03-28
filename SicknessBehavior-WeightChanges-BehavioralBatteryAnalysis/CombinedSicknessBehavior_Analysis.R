# Sickness behavior analysis across cohorts
# Jenn Kim

# load libraries
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(nlme)
library(lsmeans)

# load data
data <- read.csv("SicknessBehavior.csv") %>% select(-X)
data <- as.data.frame(data)

data$LPS.Treatment <- factor(data$LPS.Treatment, levels=c("PBS","1XLPS","2XLPS","3XLPS","4XLPS"))
data$sex <- factor(data$sex)
data$animal <- as.character(data$animal)
data$cohort <- factor(data$cohort)

data_plot <- data

# update cohort names
data_plot <- data_plot %>% mutate(cohort = 
                                    case_when(cohort=="1" ~ "RNAseq",
                                              cohort=="6" ~ "Microglia Morphology & RTqPCR",
                                              cohort=="7" ~ "Microglia Morphology & RTqPCR",
                                              cohort=="8" ~ "Behavior battery"))

### SICKNESS BEHAVIOR ###
library(ARTool)

data_plot$LPS.Treatment <- factor(data_plot$LPS.Treatment)
data_plot$cohort <- factor(data_plot$cohort)
data_plot$sex <- factor(data_plot$sex)

# model (not considering Sex since 3xLPS group is not present in females)
model <- art(Total.Score ~ LPS.Treatment, data = data_plot) 
anova = anova(model)
anova = as.data.frame(anova)
anova

write.csv(anova,"ANOVA_SicknessBehavior_ART.csv")

#posthocs 
posthoc <- art.con(model, ~LPS.Treatment, adjust="BH") %>% summary()

posthoc$Significant <- ifelse(posthoc$p.value < 0.05, "significant", "ns")

write.csv(posthoc, "SicnessBehavior_ART.Con_BH_posthocs.csv")

# boxplots for sickness behavior          
ggplot(data_plot, aes(x=LPS.Treatment, y=Total.Score, group=LPS.Treatment)) + 
  geom_boxplot(width=1, aes(group=LPS.Treatment, fill=LPS.Treatment), outlier.shape=NA) +
  geom_jitter(width = 0.2, aes(color=sex, shape=cohort), size=1.5)+
  #geom_point(position = position_dodge(width = 0.90), size=1, aes(color=sex, shape=cohort)) + 
  scale_fill_manual(values = c("lightblue", "yellow", "red","blue","purple")) +
  theme_cowplot(font_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "Condition") + #,limits = order2)+
  scale_y_continuous(name = "Total Sickness Score") +
  labs(fill="LPS Treatment", shape="Experiment", color="Sex") +
  ggtitle("Sickness behavior changes")
