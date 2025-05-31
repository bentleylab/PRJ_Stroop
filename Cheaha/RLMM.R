setwd("/Users/anaskhan/Desktop/PRJ_Stroop/results/Behavior")
df = read.csv("StroopRTdata.csv")

library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(emmeans)

df$CurrentConflict <- factor(df$CurrentConflict)
df$PreviousConflict <- factor(df$PreviousConflict)
df$CurrentConflict <- relevel(df$CurrentConflict, ref = "NoConflict")
df$PreviousConflict <- relevel(df$PreviousConflict, ref = "NoConflict")

model <- lmer(logRT ~ CurrentConflict*PreviousConflict + CurrentConflict*PC +
                (1 | Subject), data = df, REML = TRUE)

summary(model)
emm_interaction = emmeans(model, ~ CurrentConflict*PreviousConflict*PC, pbkrtest.limit = 5888)

post_hoc_tests = pairs(emm_interaction)
summary(post_hoc_tests)

emm_df <- as.data.frame(post_hoc_tests)
