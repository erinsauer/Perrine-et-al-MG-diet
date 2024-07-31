#############################################################################
############### Code for Perrine et al. Diet-MG exp     #####################
###############    developed by Erin L. Sauer           #####################
############### sections are broken down question       #####################
#############################################################################

########### load packages #############

library(tidyverse)
library(glmmTMB)
library(car)
library(lme4)
library(mgcv)

########### Experiment 1 ##############
# read in data 
E1.daily.food <- read.csv("2020 DietXDisease Daily Sumfood.csv")
E1.weekly.food <- read.csv("2020 DietXDisease Weekly Sumfood.csv")
E1.ES <- read.csv("2020 DietXDisease_Eyescore.csv")
E1.ES <- subset(E1.ES, Infection == "MG")
E1.cond <- read.csv("2020 DietXDisease_Fat and mass.csv")
E1.hem <- read.csv("2020 DietXDisease_hematocrit.csv")
E1.AB <- read.csv("2020 DietXDisease_MG Antibodies.csv")
E1.AB <- subset(E1.AB, infection == "MG")
E1.load <- read.csv("2020 DietXDisease_Pathogen load.csv")
E1.load <- subset(E1.load, Infection == "MG")
E1.WBC <- read.csv("2020 DietXDisease_WBCs.csv")
E1.Recovery <- read.csv("Exp1DaystoRecover.csv")

########### food intake analyses E1 ##############
str(E1.daily.food)
hist(E1.daily.food$dailyfood)

E1DF <- lmer(dailyfood ~ diet*infection*day + (1|ID), 
                data=E1.daily.food)
summary(E1DF)
Anova(E1DF)
E1DFS <- summary(E1DF)
str(E1DFS$coefficients)
E1DFS <- (E1DFS$coefficients)
E1DFS<- as.matrix(E1DFS)
write.csv(E1DFS, file = "E1DF_lmer.csv")
E1DFA<- as.matrix(Anova(E1DF))
colnames(E1DFA)<-c("Chi squared", "df", "p value")
write.csv(E1DFA, file = "E1DF_ANOVA.csv")

E1.daily.food$x <- paste(E1.daily.food$diet,E1.daily.food$infection)
str(E1.daily.food)
DF1<-ggplot(E1.daily.food, aes(y=dailyfood, x=day, color=x)) +
  geom_smooth()+
  #geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  xlim(-7,28)+
  ylab("Daily food intake (g)")+
  xlab("Days post inoculation")

str(E1.weekly.food)
hist(E1.weekly.food$sumfood)

E1WF <- lmer(sumfood ~ diet*infection*week + (1|ID), 
             data=E1.weekly.food)
summary(E1WF)
Anova(E1WF)
E1WFS <- summary(E1WF)
str(E1WFS$coefficients)
E1WFS <- (E1WFS$coefficients)
E1WFS<- as.matrix(E1WFS)
write.csv(E1WFS, file = "E1WF_lmer.csv")
E1WFA<- as.matrix(Anova(E1WF))
colnames(E1WFA)<-c("Chi squared", "df", "p value")
write.csv(E1WFA, file = "E1WF_ANOVA.csv")

E1.weekly.food$Treatment <- paste(E1.weekly.food$diet,E1.weekly.food$infection)
str(E1.weekly.food)
WF1 <- ggplot(E1.weekly.food, aes(y=sumfood, x=week, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="top",
        legend.text=element_text(size=15),
        legend.title=element_text(size=0))+
  xlim(0,4)+
  ylab("Weekly food intake (g)")+
  xlab("Week post inoculation")

########### eye score analysis E1 ##############
str(E1.ES)
hist(E1.ES$eyescore)

E1ES <- gamm(eyescore ~ Diet+s(Day, k=3), family=poisson,
           correlation=corCAR1(form=~Day|Bird), data = E1.ES)
plot(E1ES$gam, pages=1) #this plot lets us visually assess the fit
summary(E1ES$gam)

ES1 <- ggplot(E1.ES, aes(y=eyescore, x=Day, color=Diet)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position = "none")+
  ylab("Eye Score")+
  xlab("Days post inoculation")
########### body condition E1 ###############
str(E1.cond)
hist(E1.cond$mass)
hist(E1.cond$fat)

E1Mass <- lmer(mass ~ Diet*infection*Day + (1|Bird.ID), 
             data=E1.cond)
summary(E1Mass)
Anova(E1Mass)
E1MassA<- as.matrix(Anova(E1Mass))
colnames(E1MassA)<-c("Chi squared", "df", "p value")
write.csv(E1MassA, file = "E1Mass_ANOVA.csv")
E1MassS <- summary(E1Mass)
str(E1MassS$coefficients)
E1MassS <- (E1MassS$coefficients)
E1MassS<- as.matrix(E1MassS)
write.csv(E1MassS, file = "E1Mass_lmer.csv")

E1.cond$Treatment <- paste(E1.cond$Diet,E1.cond$infection)
str(E1.cond)
BM1 <- ggplot(E1.cond, aes(y=mass, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="right",
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))+
  #xlim(-7,28)+
  ylab("Body mass (g)")+
  xlab("Days post inoculation")

E1Fat <- glmmTMB(fat ~ Diet*infection*Day + (1|Bird.ID), 
               family=gaussian, data=E1.cond)

summary(E1Fat)
Anova(E1Fat)

E1.condMG <- subset(E1.cond, infection=="MG")
E1FatMG <- lmer(fat ~ Diet*Day + (1|Bird.ID), 
               data=E1.condMG)

summary(E1FatMG)
Anova(E1FatMG)
#emmeans(E1Fat, list(pairwise ~ Diet*infection), adjust = "tukey")

E1FatA<- as.matrix(Anova(E1Fat))
colnames(E1FatA)<-c("Chi squared", "df", "p value")
write.csv(E1FatA, file = "E1Fat_ANOVA.csv")

E1FatS <- summary(E1Fat)
str(E1FatS$coefficients)
E1FatS <- (E1FatS$coefficients)
E1FatS<- as.matrix(E1FatS)
write.csv(E1FatS, file = "E1Fat_lmer.csv")

BF1 <- ggplot(E1.cond, aes(y=fat, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(-8,28)+
  ylab("Fat score")+
  xlab("Days post inoculation")

########### antibody analysis E1 #############
str(E1.AB)
hist(E1.AB$AB)

E1AB <- glmmTMB(AB ~ diet*Day + (1|Bird.ID), 
                  family=gaussian, data=E1.AB)
summary(E1AB)
Anova(E1AB)
E1ABA<- as.matrix(Anova(E1AB))
colnames(E1ABA)<-c("Chi squared", "df", "p value")
write.csv(E1ABA, file = "E1AB_ANOVA.csv")

write.csv(summary(E1AB)$coef$cond, 
          file="E1AB_glmm.csv")

AB1 <- ggplot(E1.AB, aes(y=AB, x=Day, color=diet)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.16,.8),
        legend.text=element_text(size=15),
        legend.title=element_text(size=0))+
  #xlim(0,4)+
  ylab("MG-specific antibodies (OD)")+
  xlab("Days post inoculation")

########### MG load analysis E1 #############
str(E1.load)
hist(E1.load$log_load)

E1Load <- glmmTMB(log_load ~ Diet*day + (1|Bird), 
                 family=gaussian, data=E1.load)
summary(E1Load)
Anova(E1Load)
E1LoadA<- as.matrix(Anova(E1Load))
colnames(E1LoadA)<-c("Chi squared", "df", "p value")
write.csv(E1LoadA, file = "E1Load_ANOVA.csv")

write.csv(summary(E1Load)$coef$cond, 
          file="E1Load_glmm.csv")

PL1 <- ggplot(E1.load, aes(y=log_load, x=day, color=Diet)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position = "none")+
  xlim(2,21)+
  ylab("log10(MG load)")+
  xlab("Days post inoculation")

########### WBC analysis E1 ###############
str(E1.WBC)
hist(E1.WBC$Basophil) # very zero inflated
hist(E1.WBC$Eosinphil) # poisson 
hist(E1.WBC$Heterophil) # poisson
hist(E1.WBC$Lymphocyte) # normal
hist(E1.WBC$Monocyte) # poisson
hist(E1.WBC$H.L.Ratio) # bino

E1.WBC$Treatment <- paste(E1.WBC$Diet,E1.WBC$Infection)
str(E1.WBC)

E1Bas <- glmmTMB(Basophil ~ Diet*Infection*Day + (1|Bird.ID), 
             family=nbinom1(), data=E1.WBC)
summary(E1Bas)
Anova(E1Bas)
E1BasA<- as.matrix(Anova(E1Bas))
colnames(E1BasA)<-c("Chi squared", "df", "p value")
write.csv(E1BasA, file = "E1Bas_ANOVA.csv")

write.csv(summary(E1Bas)$coef$cond, 
          file="E1Bas_glmm.csv")

Ba1 <- ggplot(E1.WBC, aes(y=Basophil, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Basophil count")+
  xlab("Days post inoculation")

E1Eos <- glmmTMB(Eosinphil ~ Diet*Infection*Day + (1|Bird.ID), 
                 family=genpois, data=E1.WBC)
summary(E1Eos)
Anova(E1Eos)
write.csv(summary(E1Eos)$coef$cond, 
          file="E1Eos_glmm.csv")
E1Eos<- as.matrix(Anova(E1Eos))
colnames(E1Eos)<-c("Chi squared", "df", "p value")
write.csv(E1Eos, file = "E1Eos_ANOVA.csv")

Eo1 <- ggplot(E1.WBC, aes(y=Eosinphil, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Eosinphil count")+
  xlab("Days post inoculation")

E1Het <- glmmTMB(Heterophil ~ Diet*Infection*Day + (1|Bird.ID), 
                 family=genpois, data=E1.WBC)
summary(E1Het)
Anova(E1Het)
write.csv(summary(E1Het)$coef$cond, 
          file="E1Het_glmm.csv")
E1Het<- as.matrix(Anova(E1Het))
colnames(E1Het)<-c("Chi squared", "df", "p value")
write.csv(E1Het, file = "E1Het_ANOVA.csv")

He1 <- ggplot(E1.WBC, aes(y=Heterophil, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Heterophil count")+
  xlab("Days post inoculation")

E1Lym <- glmmTMB(Lymphocyte ~ Diet*Infection*Day + (1|Bird.ID), 
                 family=gaussian, data=E1.WBC)
summary(E1Lym)
Anova(E1Lym)
write.csv(summary(E1Lym)$coef$cond, 
          file="E1Lym_glmm.csv")
E1Lym<- as.matrix(Anova(E1Lym))
colnames(E1Lym)<-c("Chi squared", "df", "p value")
write.csv(E1Lym, file = "E1Lym_ANOVA.csv")

Ly1 <- ggplot(E1.WBC, aes(y=Lymphocyte, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Lymphocyte count")+
  xlab("Days post inoculation")

E1Mon <- glmmTMB(Monocyte ~ Diet*Infection*Day + (1|Bird.ID), 
                 family=genpois, data=E1.WBC)
summary(E1Mon)
Anova(E1Mon)
write.csv(summary(E1Mon)$coef$cond, 
          file="E1Mon_glmm.csv")
E1Mon<- as.matrix(Anova(E1Mon))
colnames(E1Mon)<-c("Chi squared", "df", "p value")
write.csv(E1Mon, file = "E1Mon_ANOVA.csv")

Mo1 <- ggplot(E1.WBC, aes(y=Monocyte, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Monocyte count")+
  xlab("Days post inoculation")

E1HLR <- glmmTMB(cbind(Heterophil,Lymphocyte) ~ Diet*Infection*Day + (1|Bird.ID), 
                 family=binomial, data=E1.WBC)
summary(E1HLR)
write.csv(summary(E1HLR)$coef$cond, 
          file="E1HLR_glmm.csv")
Anova(E1HLR)
E1HLR<- as.matrix(Anova(E1HLR))
colnames(E1HLR)<-c("Chi squared", "df", "p value")
write.csv(E1HLR, file = "E1HLR_ANOVA.csv")

HLR1 <- ggplot(E1.WBC, aes(y=H.L.Ratio, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.2,0.8),
        legend.text=element_text(size=15),
        legend.title=element_text(size=0))+
  #xlim(0,4)+
  ylab("Heterophil:lymphocyte ratio")+
  xlab("Days post inoculation")

########### Hematocrit E1 ##############
str(E1.hem)
E1.hem <- E1.hem %>% rename(Hem=Hematocrit....) #rename this weird column
hist(E1.hem$Hem) #normalish

E1Hem <- glmmTMB(Hem ~ diet*infection*Time.Point + (1|Bird.ID), 
                 family=gaussian, data=E1.hem)
summary(E1Hem)
Anova(E1Hem)
write.csv(summary(E1Hem)$coef$cond, 
          file="E1Hem_glmm.csv")
E1Hem<- as.matrix(Anova(E1Hem))
colnames(E1Hem)<-c("Chi squared", "df", "p value")
write.csv(E1Hem, file = "E1Hem_ANOVA.csv")

E1.hem$Treatment <- paste(E1.hem$diet,E1.hem$infection)
str(E1.hem)
Hem1 <- ggplot(E1.hem, aes(y=Hem, x=Time.Point, color=Treatment)) +
  geom_smooth(method="lm")+
  #geom_boxplot()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(-7,28)+
  ylab("Hematocrit (%)")+
  xlab("Days post inoculation")

ggplot(E1.hem, aes(y=Hem, x=Time.Point, color=diet)) +
  geom_smooth(method="lm")+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))+
  ylab("Hematocrit")
########### Days to recovery ###############
str(E1.Recovery)
hist(E1.Recovery$Days.to.recovery)
E1Reco <- glm(Days.to.recovery ~ Diet, data=E1.Recovery)
summary(E1Reco)
Anova(E1Reco)

write.csv(summary(E1Reco)$coefficients, 
          file="E1Reco_glm.csv") 
E1RecoA<- as.matrix(Anova(E1Reco))
colnames(E1RecoA)<-c("Chi squared", "df", "p value")
write.csv(E1RecoA, file = "E1Reco_ANOVA.csv")

########### Experiment 1 figures ################
#Figure 1: sum food, fat, ham, AB, ES, load
ggarrange(WF1,BF1, Hem1, AB1, ES1, PL1, ncol=3, nrow=2,
          labels=c("A","B","C","D","E","F"),
          font.label = list(size=20,face="bold"))
#saved at 1700x900

#Figure S1: daily food, mass
ggarrange(DF1,BM1, ncol=2, nrow=1,
          labels=c("A","B"),
          font.label = list(size=20,face="bold"))
#saved at 1100x400

#Figure S2: WBC
ggarrange(Ba1,Eo1,He1, Ly1,Mo1,HLR1, ncol=2, nrow=3,
          labels=c("A","B","C","D","E","F"),
          font.label = list(size=20,face="bold"))
#saved at 1100x1300

########### Experiment 2 ##############
# read in data 
E2.daily.food <- read.csv("2021 Choice_Daily food.csv")
E2.weekly.food <- read.csv("2021 Choice_Weekly Food.csv")
E2.ES <- read.csv("2021 Choice_eyescore.csv")
E2.ES <- subset(E2.ES, Infection == "MG")
E2.cond <- read.csv("2021 Choice_fat and mass.csv")
E2.hem <- read.csv("2021 Choice_hematocrit.csv")
E2.AB <- read.csv("2021 Choice_MG Antibodies.csv")
E2.AB <- subset(E2.AB, Infection == "MG")
E2.load <- read.csv("2021 Choice_pathogen load.csv")
E2.load <- subset(E2.load, Infection == "MG")
E2.WBC <- read.csv("2021 Choice_WBCs.csv")

########### food intake analyses E2 ##############
str(E2.daily.food)
hist(E2.daily.food$Lipid)
hist(E2.daily.food$Protein)
hist(E2.daily.food$Total)

E2DFL <- lmer(Lipid ~ sex*infection*Day + (1|ID), 
             data=E2.daily.food)
summary(E2DFL)
Anova(E2DFL)
E2DFAL<- as.matrix(Anova(E2DFL))
colnames(E2DFLA)<-c("Chi squared", "df", "p value")
write.csv(E2DFLA, file = "E2DFL_ANOVA.csv")

E2DFLS <- summary(E2DFL)
str(E2DFLS$coefficients)
E2DFLS <- (E2DFLS$coefficients)
E2DFLS<- as.matrix(E2DFLS)
write.csv(E2DFLS, file = "E2DFL_lmer.csv")

E2DFP <- lmer(Protein ~ sex*infection*Day + (1|ID), 
              data=E2.daily.food)
summary(E2DFP)
Anova(E2DFP)
E2DFPA<- as.matrix(Anova(E2DFP))
colnames(E2DFPA)<-c("Chi squared", "df", "p value")
write.csv(E2DFPA, file = "E2DFP_ANOVA.csv")

E2DFPS <- summary(E2DFP)
str(E2DFPS$coefficients)
E2DFPS <- (E2DFPS$coefficients)
E2DFPS<- as.matrix(E2DFPS)
write.csv(E2DFPS, file = "E2DFP_lmer.csv")

E2DFT <- lmer(Total ~ sex*infection*Day + (1|ID), 
              data=E2.daily.food)
summary(E2DFT)
Anova(E2DFT)
E2DFTA<- as.matrix(Anova(E2DFT))
colnames(E2DFTA)<-c("Chi squared", "df", "p value")
write.csv(E2DFTA, file = "E2DFT_ANOVA.csv")

E2DFTS <- summary(E2DFT)
str(E2DFTS$coefficients)
E2DFTS <- (E2DFTS$coefficients)
E2DFTS<- as.matrix(E2DFTS)
write.csv(E2DFTS, file = "E2DFT_lmer.csv")

E2.daily.food$sex <- recode_factor(E2.daily.food$sex,
                                   "M"="Male","F"="Female")
E2.daily.food$infection <- recode_factor(E2.daily.food$infection,
                                         "c"="control","mg"="MG")
E2.daily.food$Treatment <- paste(E2.daily.food$sex,E2.daily.food$infection)
str(E2.daily.food)

DF2 <- ggplot(E2.daily.food, aes(y=Total, x=Day, color=Treatment)) +
  geom_smooth()+
  #geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.2,0.85),
        legend.text=element_text(size=15),
        legend.title=element_text(size=0))+
  #xlim(0,4)+
  ylab("Total food intake (g)")+
  xlab("Days post inoculation")
DP2 <- ggplot(E2.daily.food, aes(y=Protein, x=Day, color=Treatment)) +
  geom_smooth()+
  #geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Protein diet intake (g)")+
  xlab("Days post inoculation")
DL2 <- ggplot(E2.daily.food, aes(y=Lipid, x=Day, color=Treatment)) +
  geom_smooth()+
  #geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Lipid diet intake (g)")+
  xlab("Days post inoculation")

str(E2.weekly.food)
hist(E2.weekly.food$SumLipid)
hist(E2.weekly.food$SumProtein)
hist(E2.weekly.food$SumTotal)

E2WFL <- lmer(SumLipid ~ sex*infection*Week + (1|ID), 
             data=E2.weekly.food)
summary(E2WFL)
Anova(E2WFL)
E2WFLA<- as.matrix(Anova(E2WFL))
colnames(E2WFLA)<-c("Chi squared", "df", "p value")
write.csv(E2WFLA, file = "E2WFL_ANOVA.csv")

E2WFLS <- summary(E2WFL)
str(E2WFLS$coefficients)
E2WFLS <- (E2WFLS$coefficients)
E2WFLS<- as.matrix(E2WFLS)
write.csv(E2WFLS, file = "E2WFL_lmer.csv")

E2WFP <- lmer(SumProtein ~ sex*infection*Week + (1|ID), 
              data=E2.weekly.food)
summary(E2WFP)
Anova(E2WFP)
E2WFPA<- as.matrix(Anova(E2WFP))
colnames(E2WFPA)<-c("Chi squared", "df", "p value")
write.csv(E2WFPA, file = "E2WFP_ANOVA.csv")

E2WFPS <- summary(E2WFP)
str(E2WFPS$coefficients)
E2WFPS <- (E2WFPS$coefficients)
E2WFPS<- as.matrix(E2WFPS)
write.csv(E2WFPS, file = "E2WFP_lmer.csv")

E2WFT <- lmer(SumTotal ~ sex*infection*Week + (1|ID), 
              data=E2.weekly.food)
summary(E2WFT)
Anova(E2WFT)
E2WFTA<- as.matrix(Anova(E2WFT))
colnames(E2WFTA)<-c("Chi squared", "df", "p value")
write.csv(E2WFTA, file = "E2WFT_ANOVA.csv")

E2WFTS <- summary(E2WFT)
str(E2WFTS$coefficients)
E2WFTS <- (E2WFTS$coefficients)
E2WFTS<- as.matrix(E2WFTS)
write.csv(E2WFTS, file = "E2WFT_lmer.csv")

E2.weekly.food$sex <- recode_factor(E2.weekly.food$sex,
                                   "M"="Male","F"="Female")
E2.weekly.food$infection <- recode_factor(E2.weekly.food$infection,
                                         "c"="control","mg"="MG")
E2.weekly.food$Treatment <- paste(E2.weekly.food$sex,E2.weekly.food$infection)
str(E2.weekly.food)

WF2 <- ggplot(E2.weekly.food, aes(y=SumTotal, x=Week, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.2,0.85),
        legend.text=element_text(size=15),
        legend.title=element_text(size=0))+
  #xlim(0,4)+
  ylab("Weekly total food intake (g)")+
  xlab("Weeks post inoculation")
WP2 <- ggplot(E2.weekly.food, aes(y=SumProtein, x=Week, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Weekly protein diet intake (g)")+
  xlab("Weeks post inoculation")
WL2 <- ggplot(E2.weekly.food, aes(y=SumLipid, x=Week, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Weekly lipid diet intake (g)")+
  xlab("Weeks post inoculation")

########### eye score analysis E2 ##############
str(E2.ES)
hist(E2.ES$ES)

E2ES <- gamm(ES ~ Sex+s(Day, k=3), family=poisson,
             correlation=corCAR1(form=~Day|Bird), data = E2.ES)
plot(E2ES$gam, pages=1) #this plot lets us visually assess the fit
summary(E2ES$gam)

E2.ES$Sex <- recode_factor(E2.ES$Sex,
                                   "M"="Male","F"="Female")
ES2 <- ggplot(E2.ES, aes(y=ES, x=Day, color=Sex)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position = c(0.8,0.8),
        legend.text = element_text(size=15),
        legend.title = element_text(size=0)
        )+
  ylab("Eye Score")+
  xlab("Days post inoculation")

########### body condition E2 ###############
str(E2.cond)
hist(E2.cond$mass)
hist(E2.cond$fat)

E2Mass <- lmer(mass ~ sex*Infection*Day + (1|Bird), 
               data=E2.cond)
summary(E2Mass)
Anova(E2Mass)
E2MassA<- as.matrix(Anova(E2Mass))
colnames(E2MassA)<-c("Chi squared", "df", "p value")
write.csv(E2MassA, file = "E2Mass_ANOVA.csv")

E2MassS <- summary(E2Mass)
str(E2MassS$coefficients)
E2MassS <- (E2MassS$coefficients)
E2MassS<- as.matrix(E2MassS)
write.csv(E2MassS, file = "E2Mass_lmer.csv")

E2.cond$sex <- recode_factor(E2.cond$sex,
                             "m"="Male","f"="Female")
E2.cond$Treatment <- paste(E2.cond$sex,E2.cond$Infection)
str(E2.cond)

BM2 <- ggplot(E2.cond, aes(y=mass, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="top",
        legend.text = element_text(size=15),
        legend.title = element_text(size=0))+
  #xlim(-7,28)+
  ylab("Body mass (g)")+
  xlab("Days post inoculation")

E2Fat <- lmer(fat ~ sex*Infection*Day + (1|Bird), 
              data=E2.cond)
summary(E2Fat)
Anova(E2Fat)
E2FatA<- as.matrix(Anova(E2Fat))
colnames(E2FatA)<-c("Chi squared", "df", "p value")
write.csv(E2FatA, file = "E2Fat_ANOVA.csv")

E2FatS <- summary(E2Fat)
str(E2FatS$coefficients)
E2FatS <- (E2FatS$coefficients)
E2FatS<- as.matrix(E2FatS)
write.csv(E2FatS, file = "E2Fat_lmer.csv")

BF2 <- ggplot(E2.cond, aes(y=fat, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(-7,28)+
  ylab("Fat score")+
  xlab("Days post inoculation")


########### antibody analysis E2 #############
str(E2.AB)
hist(E2.AB$AB)

E2AB <- glmmTMB(AbValue ~ Sex*Day + (1|BirdID), 
                family=gaussian, data=E2.AB)
summary(E2AB)
Anova(E2AB)
write.csv(summary(E2AB)$coef$cond, 
          file="E2AB_glmm.csv")

E2AB<- as.matrix(Anova(E2AB))
colnames(E2AB)<-c("Chi squared", "df", "p value")
write.csv(E2AB, file = "E2AB_ANOVA.csv")

AB2 <- ggplot(E2.AB, aes(y=AbValue, x=Day, color=Sex)) +
  geom_smooth(method="lm")+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.16,.9),
        legend.text=element_text(size=15),
        legend.title=element_text(size=0))+
  #xlim(0,4)+
  ylab("MG-specific antibodies (OD)")+
  xlab("Days post inoculation")

########### MG load analysis E2 #############
str(E2.load)
hist(E2.load$log.load)
E2.load <- subset(E2.load, day <= 30)

E2Load <- glmmTMB(log.load ~ sex*day + (1|Bird.ID), 
                  family=gaussian, data=E2.load)
summary(E2Load)
Anova(E2Load)
write.csv(summary(E2Load)$coef$cond, 
          file="E2Load_glmm.csv")

E2Load<- as.matrix(Anova(E2Load))
colnames(E2Load)<-c("Chi squared", "df", "p value")
write.csv(E2Load, file = "E2Load_ANOVA.csv")

PL2 <- ggplot(E2.load, aes(y=log.load, x=day, color=sex)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position = "none")+
  xlim(7,21)+
  ylab("log10(MG load)")+
  xlab("Days post inoculation")
########### WBC analysis E2 ###############
str(E2.WBC)
hist(E2.WBC$Basophil) # not analyzed, only 1 bird had >0
hist(E2.WBC$Eosinphil) # normal 
hist(E2.WBC$Heterophil) # poisson
hist(E2.WBC$Lymphocyte) # normal
hist(E2.WBC$Monocyte) # poisson
hist(E2.WBC$H.L.ratio) # bino

E2.WBC$sex <- recode_factor(E2.WBC$sex,
                                   "m"="Male","f"="Female")
E2.WBC$trt <- recode_factor(E2.WBC$trt,
                                         "C"="control")
E2.WBC$Treatment <- paste(E2.WBC$sex,E2.WBC$trt)
str(E2.WBC)

# E2Bas <- glmmTMB(Basophil ~ Diet*Infection*Day + (1|Bird.ID), 
#                  family=nbinom1(), data=E2.WBC)
# summary(E2Bas)
# Anova(E2Bas)
# E2Bas<- as.matrix(Anova(E2Bas))
# colnames(E2Bas)<-c("Chi squared", "df", "p value")
# write.csv(E2Bas, file = "E2Bas_ANOVA.csv")

E2Eos <- glmmTMB(Eosinphil ~ sex*trt*Day + (1|Bird.ID), 
                 family=gaussian, data=E2.WBC)
summary(E2Eos)
Anova(E2Eos)
write.csv(summary(E2Eos)$coef$cond, 
          file="E2Eos_glmm.csv")
E2Eos<- as.matrix(Anova(E2Eos))
colnames(E2Eos)<-c("Chi squared", "df", "p value")
write.csv(E2Eos, file = "E2Eos_ANOVA.csv")

Eo2 <- ggplot(E2.WBC, aes(y=Eosinphil, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="top",
        legend.text=element_text(size=15),
        legend.title=element_text(size=0))+
  #xlim(0,4)+
  ylab("Eosinphil count")+
  xlab("Days post inoculation")

E2Het <- glmmTMB(Heterophil ~ sex*trt*Day + (1|Bird.ID), 
                 family=genpois, data=E2.WBC)
summary(E2Het)
Anova(E2Het)
write.csv(summary(E2Het)$coef$cond, 
          file="E2Het_glmm.csv")
E2Het<- as.matrix(Anova(E2Het))
colnames(E2Het)<-c("Chi squared", "df", "p value")
write.csv(E2Het, file = "E2Het_ANOVA.csv")

He2 <- ggplot(E2.WBC, aes(y=Heterophil, x=Day, color=Treatment)) +
  geom_smooth(method="lm")+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Heterophil count")+
  xlab("Days post inoculation")

E2Lym <- glmmTMB(Lymphocyte ~ sex*trt*Day + (1|Bird.ID), 
                 family=gaussian, data=E2.WBC)
summary(E2Lym)
Anova(E2Lym)
write.csv(summary(E2Lym)$coef$cond, 
          file="E2Lym_glmm.csv")
E2Lym<- as.matrix(Anova(E2Lym))
colnames(E2Lym)<-c("Chi squared", "df", "p value")
write.csv(E2Lym, file = "E2Lym_ANOVA.csv")

Ly2 <- ggplot(E2.WBC, aes(y=Lymphocyte, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Lymphocyte count")+
  xlab("Days post inoculation")

E2Mon <- glmmTMB(Monocyte ~ sex*trt*Day + (1|Bird.ID), 
                 family=genpois, data=E2.WBC)
summary(E2Mon)
Anova(E2Mon)
write.csv(summary(E2Mon)$coef$cond, 
          file="E2Mon_glmm.csv")
E2Mon<- as.matrix(Anova(E2Mon))
colnames(E2Mon)<-c("Chi squared", "df", "p value")
write.csv(E2Mon, file = "E2Mon_ANOVA.csv")

Mo2 <- ggplot(E2.WBC, aes(y=Monocyte, x=Day, color=Treatment)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(0,4)+
  ylab("Monocyte count")+
  xlab("Days post inoculation")

E2HLR <- glmmTMB(cbind(Heterophil,Lymphocyte) ~ sex*trt*Day + (1|Bird.ID), 
                 family=binomial, data=E2.WBC)
summary(E2HLR)
Anova(E2HLR)
write.csv(summary(E2HLR)$coef$cond, 
          file="E2HLR_glmm.csv")
E2HLR<- as.matrix(Anova(E2HLR))
colnames(E2HLR)<-c("Chi squared", "df", "p value")
write.csv(E2HLR, file = "E2HLR_ANOVA.csv")

HLR2 <- ggplot(E2.WBC, aes(y=log(H.L.ratio), x=Day, color=sex)) +
  geom_smooth(method="lm")+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="top")+
  #xlim(0,4)+
  ylab("Heterophil:lymphocyte ratio")+
  xlab("Days post inoculation")


########### Hematocrit E2 ##############
str(E2.hem)
hist(E2.hem$Hematocrit) #normalish

E2Hem <- glmmTMB(Hematocrit ~ sex*Treatment*Time.Point + (1|Bird.ID), 
                 family=gaussian, data=E2.hem)
summary(E2Hem)
Anova(E2Hem)
write.csv(summary(E2Hem)$coef$cond, 
          file="E2Hem_glmm.csv")
E2Hem<- as.matrix(Anova(E2Hem))
colnames(E2Hem)<-c("Chi squared", "df", "p value")
write.csv(E2Hem, file = "E2Hem_ANOVA.csv")

E2.hem$sex <- recode_factor(E2.hem$sex,
                             "m"="Male","f"="Female")
E2.hem$Treatment2 <- paste(E2.hem$sex,E2.hem$Treatment)
str(E2.hem)
Hem2 <- ggplot(E2.hem, aes(y=Hematocrit, x=Time.Point, color=Treatment2)) +
  geom_smooth(method="lm")+
  #geom_boxplot()+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none")+
  #xlim(-7,28)+
  ylab("Hematocrit (%)")+
  xlab("Days post inoculation")


  
########### Figures experiment 2 ###########
#Figure 2: sum food - total, protein, lipid
ggarrange(WF2,WP2,WL2, ncol=3, nrow=1,
          labels=c("A","B","C"),
          font.label = list(size=20,face="bold"))
#saved at 1500x500

#Figure 3: mass, fat, ES, hem
ggarrange(BF2,BM2,Hem2,ES2, ncol=2, nrow=2,
          labels=c("A","B","C", "D"),
          font.label = list(size=20,face="bold"))
#saved at 1300x800

#Figure S3: daily - total, protein, lipid
ggarrange(DF2,DP2,DL2, ncol=3, nrow=1,
          labels=c("A","B","C"),
          font.label = list(size=20,face="bold"))
#saved at 1300x800
  
#Figure S4: AB, path load
ggarrange(AB2,PL2, ncol=2, nrow=1,
          labels=c("A","B"),
          font.label = list(size=20,face="bold"))
#saved at 700x500 
  
#Figure S5: WBC
ggarrange(Eo2,He2, Ly2,Mo2,HLR2, ncol=2, nrow=3,
          labels=c("A","B","C","D","E"),
          font.label = list(size=20,face="bold"))
#saved at 1200x1350
