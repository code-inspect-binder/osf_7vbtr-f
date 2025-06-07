######## CONFLICT TASK DATA ########
#in all analyses Age Group is coded as 0 = Adults (Reference), 1 = Early Adolescents, 2 = Mid-Adolescents, and 3 = Late Adolescents
#for Previous Trial Congruency and Current Trial Congrunecy 0 = congruent, 1 = incongruent
#for Accuracy 0 = incorrect, 1 = correct

#Code for the FIGURES in the manuscript can be found at the bottom of the file

#### Pre-process data ####
#Load the merged data files
simon <-
  read.csv(
    "simon_merge.csv",
    header = TRUE,
    sep = ",",
    na.strings = "-999",
    dec = ".",
    strip.white = TRUE
  )
flanker <-
  read.csv(
    "flanker_merge.csv",
    header = TRUE,
    sep = ",",
    na.strings = "-999",
    dec = ".",
    strip.white = TRUE
  )

#Create a variable that codes which task the data is coming from
simon <- simon[, 1:13]
simon$task <- rep(0, length(simon$subid))
flanker <- flanker[, 1:13]
flanker$task <- rep(1, length(flanker$subid))
total <- rbind(simon, flanker)

cross_task <- data.frame()

#Preprocess the data separately for each task
for (tasks in 1:2) {
  merged <- total[total$task == (tasks - 1),]
  
  #First, remove RTs below 150 ms
  
  trim1 = numeric()
  for (a in 1:length(merged$RT)) {
    if (merged$RT[a] < 150) {
      trim1[a] <- NA
    } else {
      trim1[a] <- merged$RT[a]
    }
  }
  
  merged$trim1 <- trim1
  
  #Second, standardize RTs within each participant
  merged$zrt <- rep(NA, length(merged$subid))
  for (b in 1:length(merged$RT)) {
    if (is.na(merged$trim1[b]) == TRUE) {
      merged$zrt[b] <- NA
    } else {
      merged$zrt[b] <-
        (merged$trim1[b] - mean(merged$trim1[merged$subid == merged$subid[b]], na.rm = TRUE)) / sd(merged$trim1[merged$subid == merged$subid[b]], na.rm = TRUE)
    }
  }
  
  #Third, trim based on standardized values (cutoff: 3)
  merged$trim2 <- rep(NA, length(merged$subid))
  
  for (e in 1:length(merged$trim1)) {
    if (is.na(merged$trim1[e]) == TRUE) {
      merged$trim2[e] <- NA
    } else {
      if (merged$zrt[e] < 3 && merged$zrt[e] > -3) {
        merged$trim2[e] <- merged$trim1[e]
      } else {
        merged$trim2[e] <- NA
      }
    }
  }
  
  #Fourth, re-standardize RTs after trimming
  merged$zrtf <- rep(NA, length(merged$subid))
  for (f in 1:length(merged$RT)) {
    if (is.na(merged$trim2[b]) == TRUE) {
      merged$zrtf[f] <- NA
    } else {
      merged$zrtf[f] <-
        (merged$trim2[f] - mean(merged$trim2[merged$subid == merged$subid[f]], na.rm = TRUE)) / sd(merged$trim2[merged$subid == merged$subid[f]], na.rm = TRUE)
    }
  }
  
  #Fifth, create the previous congruency and previous trial accuracy variable
  
  merged$precong <- rep(NA, length(merged$subid))
  merged$preacc <- rep(NA, length(merged$subid))
  
  for (d in 1:length(merged$trial)) {
    if (merged$trial[d] < 2) {
      merged$precong[d] <- NA
      merged$preacc[d] <- NA
    } else {
      merged$precong[d] <- merged$cong[d - 1]
      merged$preacc[d] <- merged$Acc[d - 1]
    }
  }
  
  #Finally, create a variable that tracks the number of trials across all blocks in each participant
  trials_per_block <- 97
  merged$trial_number <-
    merged$trial + (merged$block - 1) * trials_per_block
  
  trial_mean <- mean(merged$trial_num)
  merged$trial_cent <- merged$trial_num - trial_mean
  
  cross_task <- rbind(cross_task, merged)
  
}

#Add participant information 

ppt_info <-
  read.csv(
    "ppt_info.csv",
    header = TRUE,
    sep = ",",
    na.strings = "-999",
    dec = ".",
    strip.white = TRUE
  )
two_tasks <- merge(cross_task, ppt_info, by = "subid", sort = FALSE)

#Remove participants with psychiatric and neurological conditions (see readme.txt for details),
#and create data set for RT analyses by removing error and post-error trials
both_tasks <-
  subset(two_tasks, filter == 1 &
           preacc == 1 & Acc == 1 & is.na(trim2) == FALSE)

#### RT analyses in the flanker ####
flanker <- subset(both_tasks, task == 1)

library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(optimx)

flanker$subid <- as.factor(flanker$subid)
flanker$cong <- as.factor(flanker$cong)
flanker$precong <- as.factor(flanker$precong)
flanker$grs <- as.factor(flanker$grs)

#raw RT
#only terms of interest

model1 <-
  lmer(
    trim2 ~ 1 + cong * precong * grs + (1 + cong * precong | subid),
    data = flanker,
    verbose = 0,
    REML = F
  )
model2 <- lmer(
  trim2 ~ 1 + cong * precong * grs + (1 + cong | subid),
  data = flanker,
  verbose = 0,
  REML = F
)
model3 <- lmer(
  trim2 ~ 1 + cong * precong * grs + (1 | subid),
  data = flanker,
  verbose = 0,
  REML = F
)
AIC(model1, model2, model3) #model1 is preferred
Anova(model1)

#z RT
#only terms of interest

modelz1 <-
  lmer(
    zrtf ~ 1 + cong * precong * grs + (1 + cong * precong | subid),
    data = flanker,
    verbose = 0,
    REML = F,
    control = lmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb'))
  ) #model fails to converge regardless of optimizer
modelz2 <- lmer(
  zrtf ~ 1 + cong * precong * grs + (1 + cong | subid),
  data = flanker,
  verbose = 0,
  REML = F
)
modelz3 <- lmer(
  zrtf ~ 1 + cong * precong * grs + (1 | subid),
  data = flanker,
  verbose = 0,
  REML = F
)
AIC(modelz1, modelz2, modelz3) #modelz2 is preferred
Anova(modelz2)

#Post-hoc analyses of the Congruency * Age Group interaction
emms <- emmeans(modelz2, ~ cong | grs)
con <- contrast(emms, interaction = "pairwise")
pairs(con, by = NULL)

#### Equivalent ANOVAs in flanker ####
library(tidyverse)
flanker_anova <- flanker %>% 
  group_by(subid, cong, precong) %>% 
  mutate(num_gr = as.numeric(levels(grs))[grs]) %>% 
  summarise(rt = mean(trim2, na.rm = T),
            zrt = mean(zrtf, na.rm = T),
            group = mean(num_gr, na.rm = T)) %>%
  ungroup() %>% 
  mutate(subid = as.factor(subid),
         group = as.factor(group))

flanker_mod <- aov(rt ~ group*cong*precong + Error(subid/(cong*precong)) + group, data=flanker_anova)
summary(flanker_mod)

flanker_modz <- aov(zrt ~ group*cong*precong + Error(subid/(cong*precong)) + group, data=flanker_anova)
summary(flanker_modz)
  

#### RT analyses in the Simon ####
simon <- subset(both_tasks, task == 0)

library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(optimx)

simon$subid <- as.factor(simon$subid)
simon$cong <- as.factor(simon$cong)
simon$precong <- as.factor(simon$precong)
simon$grs <- as.factor(simon$grs)

#raw RT
#only terms of interest

model1 <-
  lmer(
    trim2 ~ 1 + cong * precong * grs + (1 + cong * precong | subid),
    data = simon,
    verbose = 0,
    REML = F,
    control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e5))
  )
model2 <- lmer(
  trim2 ~ 1 + cong * precong * grs + (1 + cong | subid),
  data = simon,
  verbose = 0,
  REML = F
)
model3 <- lmer(
  trim2 ~ 1 + cong * precong * grs + (1 | subid),
  data = simon,
  verbose = 0,
  REML = F
)
AIC(model1, model2, model3) #model1 is preferred
Anova(model1)

#z RT
#only terms of interest

modelz1 <-
  lmer(
    zrtf ~ 1 + cong * precong * grs + (1 + cong * precong | subid),
    data = simon,
    verbose = 0,
    REML = F,
    control = lmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb'))
  )
modelz2 <- lmer(
  zrtf ~ 1 + cong * precong * grs + (1 + cong | subid),
  data = simon,
  verbose = 0,
  REML = F
)
modelz3 <- lmer(
  zrtf ~ 1 + cong * precong * grs + (1 | subid),
  data = simon,
  verbose = 0,
  REML = F
)
AIC(modelz1, modelz2, modelz3) #modelz2 is preferred
Anova(modelz2)

#Post-hoc analyses of the Congruency * Age Group interaction
emms <- emmeans(modelz3, ~ cong | grs)
con <- contrast(emms, interaction = "pairwise")
pairs(con, by = NULL)

#### Equivalent ANOVAs in Simon ####
library(tidyverse)
simon_anova <- simon %>% 
  group_by(subid, cong, precong) %>% 
  mutate(num_gr = as.numeric(levels(grs))[grs]) %>% 
  summarise(rt = mean(trim2, na.rm = T),
            zrt = mean(zrtf, na.rm = T),
            group = mean(num_gr, na.rm = T)) %>%
  ungroup() %>% 
  mutate(subid = as.factor(subid),
         group = as.factor(group))

simon_mod <- aov(rt ~ group*cong*precong + Error(subid/(cong*precong)) + group, data=simon_anova)
summary(simon_mod)

simon_modz <- aov(zrt ~ group*cong*precong + Error(subid/(cong*precong)) + group, data=simon_anova)
summary(simon_modz)

#### RT analyses across tasks ####

library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(optimx)

both_tasks$subid <- as.factor(both_tasks$subid)
both_tasks$cong <- as.factor(both_tasks$cong)
both_tasks$precong <- as.factor(both_tasks$precong)
both_tasks$task <- as.factor(both_tasks$task)
both_tasks$grs <- as.factor(both_tasks$grs)

#raw RT
#only terms of interest

model1 <-
  lmer(
    trim2 ~ 1 + cong * precong * task * grs + (1 + cong * precong * task |
                                                 subid),
    data = both_tasks,
    verbose = 0,
    REML = F,
    control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e5))
  )
model2 <-
  lmer(
    trim2 ~ 1 + cong * precong * task * grs + (1 + cong * precong |
                                                 subid),
    data = both_tasks,
    verbose = 0,
    REML = F,
    control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e5))
  )
model3 <- lmer(
  trim2 ~ 1 + cong * precong * task * grs + (1 + cong | subid),
  data = both_tasks,
  verbose = 0,
  REML = F
)
model4 <- lmer(
  trim2 ~ 1 + cong * precong * task * grs + (1 | subid),
  data = both_tasks,
  verbose = 0,
  REML = F
)
AIC(model1, model2, model3, model4)
Anova(model1)

#z RT
#only terms of interest

modelz1 <-
  lmer(
    zrtf ~ 1 + cong * precong * task * grs + (1 + cong * precong * task |
                                                subid),
    data = both_tasks,
    verbose = 0,
    REML = F,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  )
modelz2 <-
  lmer(
    zrtf ~ 1 + cong * precong * task * grs + (1 + cong * precong |
                                                subid),
    data = both_tasks,
    verbose = 0,
    REML = F
  )
modelz3 <- lmer(
  zrtf ~ 1 + cong * precong * task * grs + (1 + cong | subid),
  data = both_tasks,
  verbose = 0,
  REML = F
)
modelz4 <- lmer(
  zrtf ~ 1 + cong * precong * task * grs + (1 | subid),
  data = both_tasks,
  verbose = 0,
  REML = F
)
AIC(modelz1, modelz2, modelz3, modelz4) #modelz1 is preferred, but it could not converge, even with different optimizers
Anova(modelz3) #or modelz1

#### Equivalent ANOVAs for both tasks ####
library(tidyverse)
both_anova <- both_tasks %>% 
  group_by(subid, cong, precong, task) %>% 
  mutate(num_gr = as.numeric(levels(grs))[grs]) %>% 
  summarise(rt = mean(trim2, na.rm = T),
            zrt = mean(zrtf, na.rm = T),
            group = mean(num_gr, na.rm = T)) %>%
  ungroup() %>% 
  mutate(subid = as.factor(subid),
         group = as.factor(group),
         task = as.factor(task))

both_mod <- aov(rt ~ group*cong*precong*task + Error(subid/(cong*precong*task)) + group, data=both_anova)
summary(both_mod)

both_modz <- aov(zrt ~ group*cong*precong*task + Error(subid/(cong*precong*task)) + group, data=both_anova)
summary(both_modz)

#### Accuracy analyses in the flanker ####

tasks_acc <- subset(two_tasks, filter == 1)
flanker_acc <- subset(tasks_acc, task == 1)

library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(optimx)

flanker_acc$subid <- as.factor(flanker_acc$subid)
flanker_acc$acc <- as.factor(flanker_acc$Acc)
flanker_acc$cong <- as.factor(flanker_acc$cong)
flanker_acc$precong <- as.factor(flanker_acc$precong)
flanker_acc$task <- as.factor(flanker_acc$task)
flanker_acc$grs <- as.factor(flanker_acc$grs)

model1 <-
  glmer(
    acc ~ 1 + cong * precong * grs + (1 + cong * precong |
                                        subid),
    data = flanker_acc,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl =
                             list(maxfun = 2e5))
  )
model2 <-
  glmer(
    acc ~ 1 + cong * precong * grs + (1 + cong |
                                        subid),
    data = flanker_acc,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun =
                                                                  2e5))
  )
model3 <-
  glmer(
    acc ~ 1 + cong * precong * grs + (1 |
                                        subid),
    data = flanker_acc,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun =
                                                                  2e5))
  )
AIC(model1, model2, model3)
Anova(model3)

#### Equivalent ANOVAs in flanker ####
library(tidyverse)
flanker_acc_anova <- flanker_acc %>%
  filter(is.na(precong) == FALSE) %>% 
  group_by(subid, cong, precong) %>% 
  mutate(num_gr = as.numeric(levels(grs))[grs],
         num_acc = as.numeric(levels(acc))[acc],) %>% 
  summarise(acc = mean(num_acc, na.rm = T),
            group = mean(num_gr, na.rm = T)) %>%
  ungroup() %>% 
  mutate(subid = as.factor(subid),
         group = as.factor(group))

flanker_acc_mod <- aov(acc ~ group*cong*precong + Error(subid/(cong*precong)) + group, data=flanker_acc_anova)
summary(flanker_acc_mod)

#### Accuracy analyses in the Simon ####

tasks_acc <- subset(two_tasks, filter == 1)
simon_acc <- subset(tasks_acc, task == 0)

library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(optimx)

simon_acc$subid <- as.factor(simon_acc$subid)
simon_acc$acc <- as.factor(simon_acc$Acc)
simon_acc$cong <- as.factor(simon_acc$cong)
simon_acc$precong <- as.factor(simon_acc$precong)
simon_acc$task <- as.factor(simon_acc$task)
simon_acc$grs <- as.factor(simon_acc$grs)

model1 <-
  glmer(
    acc ~ 1 + cong * precong * grs + (1 + cong * precong |
                                        subid),
    data = simon_acc,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl =
                             list(maxfun = 2e5))
  )
model2 <-
  glmer(
    acc ~ 1 + cong * precong * grs + (1 + cong |
                                        subid),
    data = simon_acc,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun =
                                                                  2e5))
  )
model3 <-
  glmer(
    acc ~ 1 + cong * precong * grs + (1 |
                                        subid),
    data = simon_acc,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun =
                                                                  2e5))
  )
AIC(model1, model2, model3)
Anova(model1)

#### Equivalent ANOVAs in Simon ####
library(tidyverse)
simon_acc_anova <- simon_acc %>% 
  filter(is.na(precong) == FALSE) %>% 
  group_by(subid, cong, precong) %>% 
  mutate(num_gr = as.numeric(levels(grs))[grs],
         num_acc = as.numeric(levels(acc))[acc],) %>% 
  summarise(acc = mean(num_acc, na.rm = T),
            group = mean(num_gr, na.rm = T)) %>%
  ungroup() %>% 
  mutate(subid = as.factor(subid),
         group = as.factor(group))

simon_acc_mod <- aov(acc ~ group*cong*precong + Error(subid/(cong*precong)) + group, data=simon_acc_anova)
summary(simon_acc_mod)

#### Accuracy analyses across tasks ####

tasks_acc <- subset(two_tasks, filter == 1)

library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(optimx)

tasks_acc$subid <- as.factor(tasks_acc$subid)
tasks_acc$acc <- as.factor(tasks_acc$Acc)
tasks_acc$cong <- as.factor(tasks_acc$cong)
tasks_acc$precong <- as.factor(tasks_acc$precong)
tasks_acc$task <- as.factor(tasks_acc$task)
tasks_acc$grs <- as.factor(tasks_acc$grs)

model1 <-
  glmer(
    acc ~ 1 + cong * precong * grs * task + (1 + cong * precong * task |
                                        subid),
    data = tasks_acc,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl =
                             list(maxfun = 2e5))
  )
model2 <-
  glmer(
    acc ~ 1 + cong * precong * grs * task + (1 + cong * precong |
                                        subid),
    data = tasks_acc,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun =
                                                                  2e5))
  )
model3 <-
  glmer(
    acc ~ 1 + cong * precong * grs * task + (1 + cong |
                                        subid),
    data = tasks_acc,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun =
                                                                  2e5))
  )
model4 <-
  glmer(
    acc ~ 1 + cong * precong * grs * task + (1 |
                                               subid),
    data = tasks_acc,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun =
                                                                  2e5))
  )
AIC(model1, model2, model3, model4)

#### Equivalent ANOVAs in flanker ####
library(tidyverse)
tasks_acc_anova <- tasks_acc %>%
  filter(is.na(precong) == FALSE) %>% 
  group_by(subid, cong, precong, task) %>% 
  mutate(num_gr = as.numeric(levels(grs))[grs],
         num_acc = as.numeric(levels(acc))[acc],) %>% 
  summarise(acc = mean(num_acc, na.rm = T),
            group = mean(num_gr, na.rm = T)) %>%
  ungroup() %>% 
  mutate(subid = as.factor(subid),
         group = as.factor(group),
         task = as.factor(task))

tasks_acc_mod <- aov(acc ~ group*cong*precong*task + Error(subid/(cong*precong*task)) + group, data=tasks_acc_anova)
summary(tasks_acc_mod)

######## SART DATA ########
#### Pre-process data ####
#Load the merged data file
sart <- read.csv("sart_merge.csv",header=TRUE, sep=",", na.strings="-999", dec=".", strip.white=TRUE)

#Add participant descriptives
sart <- merge(sart, ppt_info, by="subid", sort = FALSE)

sart$trialRT <- sart$RT
sart$trialRT[sart$code ==3] <- NA #this removes the RTs of probe responses (probes are coded as 3 in the Code column, 1 is Go, 2 is No Go)

#First, remove RTs below 150 ms
trim1 = numeric()
for (a in 1:length(sart$RT)) {
  if (is.na(sart$trialRT[a]) == TRUE) {
    trim1[a] <- NA
  } else if (sart$trialRT[a] < 150) {
    trim1[a] <- NA
  } else {
    trim1[a] <- sart$trialRT[a]
  }
}

sart$trim1 <- trim1 

#Second, standardize the trimmed data
sart$zrt <- rep(NA, length(sart$subid))
for (b in 1:length(sart$RT)) {
  if (is.na(sart$trim1[b]) == TRUE) {
    sart$zrt[b] <- NA
  } else {
    sart$zrt[b] <-
      (sart$trim1[b] - mean(sart$trim1[sart$subid == sart$subid[b]], na.rm = TRUE)) / sd(sart$trim1[sart$subid == sart$subid[b]], na.rm = TRUE)
  }
}


#Third, trim based on standardized values (cutoff: 3)
sart$trim2 <- numeric()

for (e in 1:length(sart$trim1)) {
  if (is.na(sart$trim1[e]) == TRUE) {
    sart$trim2[e] <- NA
  } else {
    if (sart$zrt[e] < 3 && sart$zrt[e] > -3) {
      sart$trim2[e] <- sart$trim1[e]
    } else {
      sart$trim2[e] <- NA 
    }
  }
}

#Fourth, re-standardize trimmed values
sart$zrtf = numeric()
for (f in 1:length(sart$RT)) {
  if (is.na(sart$trim2[b]) == TRUE) {
    sart$zrtf[f] <- NA
  } else {
    sart$zrtf[f] <- (sart$trim2[f] - mean(sart$trim2[sart$subid == sart$subid[f]], na.rm = TRUE)) / sd(sart$trim2[sart$subid == sart$subid[f]], na.rm = TRUE) 
  }
}

#Finally, create trial number variable
trials_per_block <- 131
sart$trial_number <- sart$trial + (sart$block-1)*trials_per_block 

trial_mean <- mean(sart$trial_num)
sart$trial_cent <- sart$trial_num - trial_mean

#Remove participants with psychiatric/neurological conditions, plus an additional ppt,
#Subject 107 because the fire alarm went off during his SART
sart <- subset(sart, filter == 1 & subid != 107)

#### RT analyses in the SART ####
sart$subid <- as.factor(sart$subid)
sart$grs <- as.factor(sart$grs)
sart$sex <- as.factor(sart$sex)

#Create a binary trial type variable
#0 = Go, 1 = No Go
sart$trial_type <- sart$code
sart$trial_type[sart$code == 3] <- NA
sart$trial_type <- as.factor(sart$trial_type-1)


#raw RT
model1 <- lmer(trim2 ~ 1 + trial_type*grs + (1|subid), 
               data = sart, verbose = 0, REML = F)
model2 <- lmer(trim2 ~ 1 + trial_type*grs + (1+trial_type|subid), 
               data = sart, verbose = 0, REML = F)
AIC(model1, model2)
Anova(model2)

#zRT
modelz1 <- lmer(zrtf ~ 1 + trial_type*grs + (1|subid), 
               data = sart, verbose = 0, REML = F)
modelz2 <- lmer(zrtf ~ 1 + trial_type*grs + (1+trial_type|subid), 
               data = sart, verbose = 0, REML = F)
AIC(modelz1, modelz2)
Anova(modelz1)

#### Equivalent ANOVAs in the SART ####
#Some participants made no errors, they need to be dropped 
#so the ANOVA runs smoothly
library(tidyverse)
sart_filter <- sart %>% 
  filter(trial_type == 1) %>%
  group_by(subid) %>% 
  summarise(nogo_acc = mean(Acc, na.rm = T)) %>% 
  mutate(anova_filt = case_when(nogo_acc == 1 ~ 0,
                                TRUE ~ 1)) %>% 
  select(subid, anova_filt)

sart_anova <- sart %>% 
  filter(is.na(trial_type) == FALSE) %>% 
  group_by(subid, trial_type) %>% 
  mutate(num_gr = as.numeric(levels(grs))[grs]) %>% 
  summarise(rt = mean(trim2, na.rm = T),
            zrt = mean(zrtf, na.rm = T),
            group = mean(num_gr, na.rm = T)) %>%
  ungroup() %>% 
    inner_join(., sart_filter, by = "subid") %>%
  filter(anova_filt == 1) %>% 
  mutate(subid = as.factor(subid),
         group = as.factor(group))

sart_mod <- aov(rt ~ group*trial_type + Error(subid/(trial_type)) + group, data=sart_anova)
summary(sart_mod)

sart_modz <- aov(zrt ~ group*trial_type + Error(subid/(trial_type)) + group, data=sart_anova)
summary(sart_modz)

#### Accuracy analyses in the SART ####
sart$grs <- as.factor(sart$grs)
sart$sex <- as.factor(sart$sex)

#Create a binary trial type variable
sart$trial_type <- sart$code
sart$trial_type[sart$code == 3] <- NA
sart$trial_type <- as.factor(sart$trial_type-1)

#Create the accuracy variable
sart$acc <- sart$Acc
sart$acc[sart$code == 3] <- NA #remove the accuracy of probe responses (coded as 999)
sart$acc <- as.factor(sart$acc)

model1 <-
  glmer(
    acc ~ 1 + trial_type * grs + (1 |
                                    subid),
    data = sart,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun =
                                                                  2e5))
  )
model2 <-
  glmer(
    acc ~ 1 + trial_type * grs + (1 + trial_type |
                                    subid),
    data = sart,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun =
                                                                  2e5))
  )
AIC(model1, model2)
Anova(model2)

#Calculate odds ratios and confidence intervals for the ORs
OR <- exp(fixef(model2))
CI <- exp(confint(model2, parm = "trial_type1"))

#Alternative way of calculating CIs, as the profile way (above) fails for certain analyses
exp(summary(model2)$coefficients["trial_type1", 1] + qnorm(c(0.025, 0.5, 0.975)) * summary(model2)$coefficients["trial_type1", 2])

#Post-hoc analyses
#Group differences
emms1 <- emmeans(model2, ~ grs)
pairs(emms1, by = NULL)
#Trial Type * Age Group interaction
emms2 <- emmeans(model2, ~ trial_type | grs)
con1 <- contrast(emms2, interaction = "pairwise")
contrast(emms2, interaction = "pairwise", by = NULL)
pairs(con1, by = NULL)

#### Equivalent ANOVAs in the SART ####
library(tidyverse)
sart_acc_anova <- sart %>% 
  filter(is.na(trial_type) == FALSE) %>% 
  group_by(subid, trial_type) %>% 
  mutate(num_gr = as.numeric(levels(grs))[grs],
         num_acc = as.numeric(levels(acc))[acc],) %>% 
  summarise(acc = mean(num_acc, na.rm = T),
            group = mean(num_gr, na.rm = T)) %>%
  ungroup() %>% 
  mutate(subid = as.factor(subid),
         group = as.factor(group))

sart_acc_mod <- aov(acc ~ group*trial_type + Error(subid/(trial_type)) + group, data=sart_acc_anova)
summary(sart_acc_mod)

#### Mind-wandering & Age ####
#Create four binary variables for each thought report category
sart$p1 <- rep(0, length(sart$resp))
sart$p2 <- rep(0, length(sart$resp))
sart$p3 <- rep(0, length(sart$resp))
sart$p4 <- rep(0, length(sart$resp))

#Create a fifth category that codes both MW (with/without awareness)
#categories against all other categories
sart$mw <- rep(0, length(sart$resp))

for (f in 1:length(sart$subid)) {
  if (sart$resp[f] == 97) {
    sart$p1[f] <- 1
    sart$mw[f] <- 0
  } else if (sart$resp[f] == 98) {
    sart$p2[f] <- 1
    sart$mw[f] <- 0
  } else if (sart$resp[f] == 99) {
    sart$p3[f] <- 1
    sart$mw[f] <- 1
  } else if (sart$resp[f] == 100) {
    sart$p4[f] <- 1
    sart$mw[f] <- 1
  }
}

#Create the data set for the binary logistic regressions
#Only keep probe trials (this will result in 10 rows per participant)
probes <- subset(sart, code == 3)

#State level (SART probes)

probes$grs <- as.factor(probes$grs)

mod1 <- glmer(p1 ~ 1 + grs + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p2 ~ 1 + grs + (1|subid), data = probes, family = binomial)
mod3 <- glmer(p3 ~ 1 + grs + (1|subid), data = probes, family = binomial)
mod4 <- glmer(p4 ~ 1 + grs + (1|subid), data = probes, family = binomial)
mod5 <- glmer(mw ~ 1 + grs + (1|subid), data = probes, family = binomial)

Anova(mod1)
Anova(mod2)
Anova(mod3)
Anova(mod4)
Anova(mod5)

emms <- emmeans(mod5, "grs")
pairs(emms)

confint(pairs(emms),type = "response")

#### Equivalent ANOVA for age effect on MW ####
library(tidyverse)
mw_age_anova <- probes %>% 
  group_by(subid) %>% 
  mutate(num_gr = as.numeric(levels(grs))[grs]) %>% 
  summarise(on = mean(p1, na.rm = T),
            space = mean(p2, na.rm = T),
            zone = mean(p3, na.rm = T),
            tune = mean(p4, na.rm = T),
            mw = mean(mw, na.rm = T),
            group = mean(num_gr, na.rm = T)) %>%
  ungroup() %>% 
  mutate(subid = as.factor(subid),
         group = as.factor(group))

mw_age_mod <- aov(mw ~ group, data=mw_age_anova)
summary(mw_age_mod)


#### Table 3A - MW and Age ####
tab <- c(1:(5*2))
dim(tab) <- c(5,2)
row.names(tab) <- c("On-task", "Space Out", "Zone Out", "Tune Out", "Overall MW")
colnames(tab) <- c("Chisq", "p")
tab["On-task", "Chisq"] <- round(Anova(mod1)$Chisq, 2)
tab["On-task", "p"] <- round(Anova(mod1)$"Pr(>Chisq)", 3)
tab["Space Out", "Chisq"] <- round(Anova(mod2)$Chisq, 2)
tab["Space Out", "p"] <- round(Anova(mod2)$"Pr(>Chisq)", 3)
tab["Zone Out", "Chisq"] <- round(Anova(mod3)$Chisq, 2)
tab["Zone Out", "p"] <- round(Anova(mod3)$"Pr(>Chisq)", 3)
tab["Tune Out", "Chisq"] <- round(Anova(mod4)$Chisq, 2)
tab["Tune Out", "p"] <- round(Anova(mod4)$"Pr(>Chisq)", 3)
tab["Overall MW", "Chisq"] <- round(Anova(mod5)$Chisq, 2)
tab["Overall MW", "p"] <- round(Anova(mod5)$"Pr(>Chisq)", 3)


#Trait level (MWQ)

library(tidyr)
library(dplyr)

sart$grs <- as.numeric(levels(sart$grs))[sart$grs]
traitmw <- sart %>%
  group_by(subid) %>%
  summarise(
    mwDelib = mean(mwDelib, na.rm = TRUE),
    mwSpont = mean(mwSpont, na.rm = TRUE),
    grs = mean(grs, na.rm = TRUE)
  )
traitmw <-
  gather(traitmw, mwtype, score, mwDelib:mwSpont, factor_key = TRUE)
traitmw$grs <- as.factor(traitmw$grs)
contrasts(traitmw$grs) <- contr.poly(4)

traitmw_anova <-
  aov(score ~ grs * mwtype + Error(subid / mwtype), data = traitmw)
summary(traitmw_anova)

#### Alternative indices of MW (Go accuracy, No Go accuracy, Go RT variability) ####
sart$gort <- sart$trim2
sart$gort[sart$code != 1] <- NA
sart$nogort <- sart$trim2
sart$nogort[sart$code != 2] <- NA
sart$goacc <- sart$Acc
sart$goacc[sart$code != 1] <- NA
sart$nogoacc <- sart$Acc
sart$nogoacc[sart$code != 2] <- NA

sartsub <- sart %>%
  group_by(subid) %>%
  summarise(
    gomean = mean(gort, na.rm = TRUE),
    gosd = sd(gort, na.rm = TRUE),
    nogomean = mean(nogort, na.rm = TRUE),
    goacc_sub = mean(goacc, na.rm = TRUE),
    nogoacc_sub = mean(nogoacc, na.rm = TRUE),
    mwDelib = mean(mwDelib, na.rm = TRUE),
    mwSpont = mean(mwSpont, na.rm = TRUE)
  )

##Calculate Go RT coefficient of variation (CV)
##to index intraindividual variability
sartsub$go_cv <- sartsub$gosd / sartsub$gomean
sart <- merge(sart, sartsub[,c(1,5,6,9)], by = "subid")

##Are they related to MW?
probes <- subset(sart, code == 3)

#Are behavioural indices related to self-reports?
#check significance after each block of models using the list of ANOVAs below
mod1 <- glmer(p1 ~ 1 + goacc_sub + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p2 ~ 1 + goacc_sub + (1|subid), data = probes, family = binomial)
mod3 <- glmer(p3 ~ 1 + goacc_sub + (1|subid), data = probes, family = binomial)
mod4 <- glmer(p4 ~ 1 + goacc_sub + (1|subid), data = probes, family = binomial)
mod5 <- glmer(mw ~ 1 + goacc_sub + (1|subid), data = probes, family = binomial)

mod1 <- glmer(p1 ~ 1 + nogoacc_sub + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p2 ~ 1 + nogoacc_sub + (1|subid), data = probes, family = binomial)
mod3 <- glmer(p3 ~ 1 + nogoacc_sub + (1|subid), data = probes, family = binomial)
mod4 <- glmer(p4 ~ 1 + nogoacc_sub + (1|subid), data = probes, family = binomial)
mod5 <- glmer(mw ~ 1 + nogoacc_sub + (1|subid), data = probes, family = binomial)

mod1 <- glmer(p1 ~ 1 + go_cv + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p2 ~ 1 + go_cv + (1|subid), data = probes, family = binomial)
mod3 <- glmer(p3 ~ 1 + go_cv + (1|subid), data = probes, family = binomial)
mod4 <- glmer(p4 ~ 1 + go_cv + (1|subid), data = probes, family = binomial)
mod5 <- glmer(mw ~ 1 + go_cv + (1|subid), data = probes, family = binomial)

Anova(mod1)
Anova(mod2)
Anova(mod3)
Anova(mod4)
Anova(mod5)

#Is trait MW related to state MW?
mod1 <- glmer(p1 ~ 1 + mwDelib + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p2 ~ 1 + mwDelib + (1|subid), data = probes, family = binomial)
mod3 <- glmer(p3 ~ 1 + mwDelib + (1|subid), data = probes, family = binomial)
mod4 <- glmer(p4 ~ 1 + mwDelib + (1|subid), data = probes, family = binomial)
mod5 <- glmer(mw ~ 1 + mwDelib + (1|subid), data = probes, family = binomial)

mod1 <- glmer(p1 ~ 1 + mwSpont + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p2 ~ 1 + mwSpont + (1|subid), data = probes, family = binomial)
mod3 <- glmer(p3 ~ 1 + mwSpont + (1|subid), data = probes, family = binomial)
mod4 <- glmer(p4 ~ 1 + mwSpont + (1|subid), data = probes, family = binomial)
mod5 <- glmer(mw ~ 1 + mwSpont + (1|subid), data = probes, family = binomial)

Anova(mod1)
Anova(mod2)
Anova(mod3)
Anova(mod4)
Anova(mod5)

OR1 <- exp(fixef(mod1))
CI1 <- exp(confint(mod1))
OR4 <- exp(fixef(mod4))
CI4 <- exp(confint(mod4))

#Are behavioural indices related to trait MW?
#Variability
cor.test(sartsub$go_cv, sartsub$mwDelib)
cor.test(sartsub$go_cv, sartsub$mwSpont)

#Accuracy
mod_delib <-
  glmer(
    acc ~ 1 + mwDelib + (1 |
                           subid),
    data = sart,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
mod_spont <-
  glmer(
    acc ~ 1 + mwSpont + (1 |
                           subid),
    data = sart,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )


Anova(mod_delib)
Anova(mod_spont)

##Does RT variability (CV) change as a function of age?
#if the Age Group variable has already been converted to a factor before
#reconvert to numeric for averaging across trials
sart$grs <- as.numeric(levels(sart$grs))[sart$grs]

cvsub <- sart %>%
  group_by(subid) %>%
  summarise(
    grs = mean(grs, na.rm = TRUE),
    cv = mean(go_cv, na.rm = TRUE)
  )
cvsub$grs <- as.factor(cvsub$grs)
contrasts(cvsub$grs)<-contr.poly(4)
cvanova <- aov(cv ~ grs, data = cvsub)
summary(cvanova)
posthoc <- TukeyHSD(cvanova, 'grs', conf.level=0.95)


#### Mind-wandering & Cognitive Control ####
library(tidyverse)
#Change here which conflict task you get the CSE from
conflict_task <- subset(both_tasks, task == 1) #0 is Simon, 1 is flanker

#Get both raw and standardized CSEs per participant
csedat <- conflict_task %>%
  group_by(subid, precong, cong) %>%
  summarise(mean1 = mean(trim2, na.rm = TRUE))

cse_wide <- csedat %>%
  unite(conds, c("precong","cong")) %>%
  spread(conds,mean1)

cse_wide$cse <- (cse_wide$'0_1' - cse_wide$'0_0') - (cse_wide$'1_1' - cse_wide$'1_0')
cse_wide$ce <- ((cse_wide$'0_1' + cse_wide$'1_1')/2) - ((cse_wide$'0_0' + cse_wide$'1_0')/2)

csedat <- conflict_task %>%
  group_by(subid, precong, cong) %>%
  summarise(meanz = mean(zrtf, na.rm = TRUE))

cse_widez <- csedat %>%
  unite(conds, c("precong","cong")) %>%
  spread(conds,meanz)

cse_widez$csez <- (cse_widez$'0_1' - cse_widez$'0_0') - (cse_widez$'1_1' - cse_widez$'1_0')
cse_widez$cez <- ((cse_widez$'0_1' + cse_widez$'1_1')/2) - ((cse_widez$'0_0' + cse_widez$'1_0')/2)

cse_wide <- merge(cse_wide, cse_widez[,c("subid","csez","cez")], by = "subid")

#Add CSE variables to the MW probe data set
probes <- merge(probes, cse_wide[,c("subid","cse","ce","csez","cez")], by = "subid")

#Analyses - CSE
mod1 <- glmer(p1 ~ 1 + cse + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p2 ~ 1 + cse + (1|subid), data = probes, family = binomial)
mod3 <- glmer(p3 ~ 1 + cse + (1|subid), data = probes, family = binomial)
mod4 <- glmer(p4 ~ 1 + cse + (1|subid), data = probes, family = binomial)
mod5 <- glmer(mw ~ 1 + cse + (1|subid), data = probes, family = binomial)

Anova(mod1)
Anova(mod2)
Anova(mod3)
Anova(mod4)
Anova(mod5)

modz1 <- glmer(p1 ~ 1 + csez + (1|subid), data = probes, family = binomial)
modz2 <- glmer(p2 ~ 1 + csez + (1|subid), data = probes, family = binomial)
modz3 <- glmer(p3 ~ 1 + csez + (1|subid), data = probes, family = binomial)
modz4 <- glmer(p4 ~ 1 + csez + (1|subid), data = probes, family = binomial)
modz5 <- glmer(mw ~ 1 + csez + (1|subid), data = probes, family = binomial)

Anova(modz1)
Anova(modz2)
Anova(modz3)
Anova(modz4)
Anova(modz5)

#### Table 3B - MW and CSE - raw ####
tab <- c(1:(5*2))
dim(tab) <- c(5,2)
row.names(tab) <- c("On-task", "Space Out", "Zone Out", "Tune Out", "Overall MW")
colnames(tab) <- c("Chisq", "p")
tab["On-task", "Chisq"] <- round(Anova(mod1)$Chisq, 2)
tab["On-task", "p"] <- round(Anova(mod1)$"Pr(>Chisq)", 3)
tab["Space Out", "Chisq"] <- round(Anova(mod2)$Chisq, 2)
tab["Space Out", "p"] <- round(Anova(mod2)$"Pr(>Chisq)", 3)
tab["Zone Out", "Chisq"] <- round(Anova(mod3)$Chisq, 2)
tab["Zone Out", "p"] <- round(Anova(mod3)$"Pr(>Chisq)", 3)
tab["Tune Out", "Chisq"] <- round(Anova(mod4)$Chisq, 2)
tab["Tune Out", "p"] <- round(Anova(mod4)$"Pr(>Chisq)", 3)
tab["Overall MW", "Chisq"] <- round(Anova(mod5)$Chisq, 2)
tab["Overall MW", "p"] <- round(Anova(mod5)$"Pr(>Chisq)", 3)

#### Table 3C - MW and CSE - zRT ####
tab <- c(1:(5*2))
dim(tab) <- c(5,2)
row.names(tab) <- c("On-task", "Space Out", "Zone Out", "Tune Out", "Overall MW")
colnames(tab) <- c("Chisq", "p")
tab["On-task", "Chisq"] <- round(Anova(modz1)$Chisq, 2)
tab["On-task", "p"] <- round(Anova(modz1)$"Pr(>Chisq)", 3)
tab["Space Out", "Chisq"] <- round(Anova(modz2)$Chisq, 2)
tab["Space Out", "p"] <- round(Anova(modz2)$"Pr(>Chisq)", 3)
tab["Zone Out", "Chisq"] <- round(Anova(modz3)$Chisq, 2)
tab["Zone Out", "p"] <- round(Anova(modz3)$"Pr(>Chisq)", 3)
tab["Tune Out", "Chisq"] <- round(Anova(modz4)$Chisq, 2)
tab["Tune Out", "p"] <- round(Anova(modz4)$"Pr(>Chisq)", 3)
tab["Overall MW", "Chisq"] <- round(Anova(modz5)$Chisq, 2)
tab["Overall MW", "p"] <- round(Anova(modz5)$"Pr(>Chisq)", 3)

#Analyses - Congruency effect

#NOTE: for this first set of analyses, put scale() around the "ce" variable
#for the models to converge. This is ONLY necessary for the SIMON task.
mod1 <- glmer(p1 ~ 1 + scale(ce) + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p2 ~ 1 + scale(ce) + (1|subid), data = probes, family = binomial)
mod3 <- glmer(p3 ~ 1 + scale(ce) + (1|subid), data = probes, family = binomial)
mod4 <- glmer(p4 ~ 1 + scale(ce) + (1|subid), data = probes, family = binomial)
mod5 <- glmer(mw ~ 1 + scale(ce) + (1|subid), data = probes, family = binomial)

Anova(mod1)
Anova(mod2)
Anova(mod3)
Anova(mod4)
Anova(mod5)

modz1 <- glmer(p1 ~ 1 + cez + (1|subid), data = probes, family = binomial)
modz2 <- glmer(p2 ~ 1 + cez + (1|subid), data = probes, family = binomial)
modz3 <- glmer(p3 ~ 1 + cez + (1|subid), data = probes, family = binomial)
modz4 <- glmer(p4 ~ 1 + cez + (1|subid), data = probes, family = binomial)
modz5 <- glmer(mw ~ 1 + cez + (1|subid), data = probes, family = binomial)

Anova(modz1)
Anova(modz2)
Anova(modz3)
Anova(modz4)
Anova(modz5)


#### Alternative indices of MW ####
sartsub <- merge(sartsub, cse_wide[,c(1,6,7)], by = "subid")
#Trait MW
cor.test(sartsub$cse, sartsub$mwDelib, method = "pearson")
cor.test(sartsub$csez, sartsub$mwDelib, method = "pearson")

cor.test(sartsub$cse, sartsub$mwSpont, method = "pearson")
cor.test(sartsub$csez, sartsub$mwSpont, method = "pearson")

#Variability
cor.test(sartsub$cse, sartsub$go_cv, method = "pearson")
cor.test(sartsub$csez, sartsub$go_cv, method = "pearson")

#Accuracy
sart <- merge(sart, cse_wide[,c(1,6,7)], by = "subid")

mod <-
  glmer(
    acc ~ 1 + cse + (1 |
                           subid),
    data = sart,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
modz <-
  glmer(
    acc ~ 1 + csez + (1 |
                           subid),
    data = sart,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )

Anova(mod)
Anova(modz)

#### SUPPLEMENTARY MODELS ####
#### Controlling for potential confounds ####
#Confounds were controlled for using the following logic:

#Set which task you want to double-check here:
task <- flanker
#Next, define the control variable
#Potential control variables include:
#Sex, IQ, Positive or Negative Mood
#the corresponding variables are: sex, IQ, panas_pos1/2, panas_neg1/2
#(panas_pos1 and panas_neg1 refer to the results of the PANAS completed BEFORE the start of the conflict tasks,
#panas_pos2 and panas_neg2 refer to the results of the PANAS completed AFTER the flanker and the Simon, but BEFORE
#the SART, thus we used these latter variables to control for mood in the SART)
task$control_var <- task$sex
#Alternatively, to control for fatigue and practice effects (see Aschenbrenner & Balota, 2017), we added a
#Trial Number by Age Group interaction (trial_cent*grs) term to the models
#In case the control variable needs to be a factor (e.g., Sex):
task$control_var <- as.factor(task$control_var)

#Note: in SART analyses the cong and precong variables have to be replaced by trial_type
model1 <- lmer(trim2 ~ 1 + cong*precong*grs + control_var + (1+cong*precong|subid), 
               data = task, verbose = 0, REML = F)
model2 <- lmer(trim2 ~ 1 + cong*precong*grs + control_var + (1+cong|subid), 
               data = task, verbose = 0, REML = F)
model3 <- lmer(trim2 ~ 1 + cong*precong*grs*control_var + (1+cong*precong|subid), 
               data = task, verbose = 0, REML = F)
model4 <- lmer(trim2 ~ 1 + cong*precong*grs*control_var + (1+cong|subid), 
               data = task, verbose = 0, REML = F)
AIC(model1, model2, model3, model4)
Anova(model1)

#These analyses can be repeated using standardized RT as an outcome (by changing "trim2" to "zrtf"),
#and using age as a continuous variable as opposed to categorical (by changing "grs" to "age_c")

#Finally, controlling potential confounds in the most important MW-related findings:
#Age effect:
probes$control_var <- probes$sex
probes$control_var <- as.factor(probes$control_var)
mod0 <- glmer(p4 ~ 1 + grs + (1|subid), data = probes, family = binomial)
mod1 <- glmer(p4 ~ 1 + grs + control_var + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p4 ~ 1 + grs*control_var + (1|subid), data = probes, family = binomial)

AIC(mod0, mod1, mod2)

#CSE-MW relationship:
#Age Group and Trial Number could be especially important here
probes$control_var <- probes$trial_cent
probes$control_var <- as.factor(probes$control_var)

#Trial Number may need to be rescaled, using the scale() function
probes$control_var <- scale(probes$control_var)

mod0 <- glmer(p1 ~ 1 + cse + (1|subid), data = probes, family = binomial)
mod1 <- glmer(p1 ~ 1 + cse + control_var + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p1 ~ 1 + cse*control_var + (1|subid), data = probes, family = binomial)

AIC(mod0, mod1, mod2)

mod0 <- glmer(p1 ~ 1 + csez + (1|subid), data = probes, family = binomial)
mod1 <- glmer(p1 ~ 1 + csez + control_var + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p1 ~ 1 + csez*control_var + (1|subid), data = probes, family = binomial)

AIC(mod0, mod1, mod2)

mod0 <- glmer(p4 ~ 1 + cse + (1|subid), data = probes, family = binomial)
mod1 <- glmer(p4 ~ 1 + cse + control_var + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p4 ~ 1 + cse*control_var + (1|subid), data = probes, family = binomial)

AIC(mod0, mod1, mod2)

mod0 <- glmer(p4 ~ 1 + csez + (1|subid), data = probes, family = binomial)
mod1 <- glmer(p4 ~ 1 + csez + control_var + (1|subid), data = probes, family = binomial)
mod2 <- glmer(p4 ~ 1 + csez*control_var + (1|subid), data = probes, family = binomial)

AIC(mod0, mod1, mod2)

#### Figures ####
library(ggplot2)
#### Fig. 1 ####
##The effect of congruency as a function of previous trial congruency, or the congruency sequence effect in raw RT in the flanker task
th <- theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            strip.background =element_rect(fill="white", colour = "black"),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.text = element_text(color = "black", size = 8),
            legend.position = "right",
            legend.title = element_text(size=10),
            legend.text = element_text(size=8),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            text = element_text(face = "bold", size = 12))


#Convert group variable to numeric so it will be easier to group it
flanker$num_grs <- as.numeric(levels(flanker$grs))[flanker$grs]

interim_data <- flanker %>%
  mutate(prevcong = case_when(precong == 1L ~ "Incongruent",
                              precong == 0L ~ "Congruent"),
         congruent = case_when(cong == 1L ~ "Incongruent",
                               cong == 0L ~ "Congruent")) %>%
  group_by(subid, congruent, prevcong) %>%
  summarize(mean_rt = mean(trim2, na.rm = TRUE),
            group = mean(num_grs, na.rm = TRUE))

interim_data$group <- factor(interim_data$group,
                             levels = c(1,2,3,0),
                             labels = c("Early adolescents", "Mid-adolescents", "Late adolescents", "Adults")) 

CSE_plot_data <- interim_data %>%
  group_by(group, congruent, prevcong) %>%
  summarize(N = n(),
            mean = mean(mean_rt, na.rm = TRUE),
            sd = sd(mean_rt, na.rm = TRUE),
            se = sd / sqrt(N))


CSE_plot <- ggplot(CSE_plot_data, aes(x =prevcong, y = mean,
                                      group = congruent, shape = congruent)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_line() +
  geom_point(size = 2) +
  scale_shape_manual(values=c(4, 16)) +
  xlab("Congruency of Trial N-1")+
  ylab("Reaction time (ms)") +
  guides(shape = guide_legend(title="Congruency of \n Trial N")) +
  th +
  theme(legend.title.align=0.5) +
  theme(legend.position = c(0.9, 0.84)) +
  facet_wrap(~ group, ncol = 4)

CSE_plot

#### Fig. 2 ####
##The effect of congruency as a function of previous trial congruency, or the congruency sequence effect in raw RT in the Simon task
th <- theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            strip.background =element_rect(fill="white", colour = "black"),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.text = element_text(color = "black", size = 8),
            legend.position = "right",
            legend.title = element_text(size=10),
            legend.text = element_text(size=8),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            text = element_text(face = "bold", size = 12))


#Convert group variable to numeric so it will be easier to group it
simon$num_grs <- as.numeric(levels(simon$grs))[simon$grs]

interim_data <- simon %>%
  mutate(prevcong = case_when(precong == 1L ~ "Incongruent",
                              precong == 0L ~ "Congruent"),
         congruent = case_when(cong == 1L ~ "Incongruent",
                               cong == 0L ~ "Congruent")) %>%
  group_by(subid, congruent, prevcong) %>%
  summarize(mean_rt = mean(trim2, na.rm = TRUE),
            group = mean(num_grs, na.rm = TRUE))

interim_data$group <- factor(interim_data$group,
                             levels = c(1,2,3,0),
                             labels = c("Early adolescents", "Mid-adolescents", "Late adolescents", "Adults")) 

CSE_plot_data <- interim_data %>%
  group_by(group, congruent, prevcong) %>%
  summarize(N = n(),
            mean = mean(mean_rt, na.rm = TRUE),
            sd = sd(mean_rt, na.rm = TRUE),
            se = sd / sqrt(N))


CSE_plot <- ggplot(CSE_plot_data, aes(x =prevcong, y = mean,
                                      group = congruent, shape = congruent)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_line() +
  geom_point(size = 2) +
  scale_shape_manual(values=c(4, 16)) +
  xlab("Congruency of Trial N-1")+
  ylab("Reaction time (ms)") +
  guides(shape = guide_legend(title="Congruency of \n Trial N")) +
  th +
  theme(legend.title.align=0.5) +
  theme(legend.position = c(0.9, 0.84)) +
  facet_wrap(~ group, ncol = 4)

CSE_plot

#### Fig. 3 ####
##Box-plots of the frequencies of different categories of thought reports across age groups during the SART
if (is.factor(probes$grs)) {
  probes$num_grs <- as.numeric(levels(probes$grs))[probes$grs]
} else {
  probes$num_grs <- probes$grs
}
subprobe_data <- probes %>%
  complete(subid) %>%
  group_by(subid) %>% 
  summarize(p1 = sum(p1, na.rm = TRUE),
            p2 = sum(p2, na.rm = TRUE),
            p3 = sum(p3, na.rm = TRUE),
            p4 = sum(p4, na.rm = TRUE),
            grs = mean(num_grs, na.rm = TRUE))
subprobe_data$grs <- factor(subprobe_data$grs,
                             levels = c(1,2,3,0),
                             labels = c("Early adolescents", "Mid-adolescents", "Late adolescents", "Adults"))

plot1 <- ggplot(subprobe_data, aes(x=grs, y=p1)) + 
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(face = "bold", size = 12)) +
  theme(axis.text = element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  labs(x = "Age Groups", y = "On-Task Report Frequency") +
  ylim(0, 10)
plot2 <- ggplot(subprobe_data, aes(x=grs, y=p2)) + 
  geom_boxplot() +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(text = element_text(face = "bold", size = 12)) +
  theme(axis.text = element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  labs(y = "Space Out Frequency") +
  ylim(0, 10)
plot3 <- ggplot(subprobe_data, aes(x=grs, y=p3)) + 
  geom_boxplot() +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(text = element_text(face = "bold", size = 12)) +
  theme(axis.text = element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  labs(y = "Zone Out Frequency") +
  ylim(0, 10)
plot4 <- ggplot(subprobe_data, aes(x=grs, y=p4)) + 
  geom_boxplot() +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(text = element_text(face = "bold", size = 12)) +
  theme(axis.text = element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  labs(y = "Tune Out Frequency") +
  ylim(0, 10)

library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot4), ggplotGrob(plot3), ggplotGrob(plot2), ggplotGrob(plot1), size = "last"))