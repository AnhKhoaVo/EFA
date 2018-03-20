#Raw and Scaled scores for MS52 and MS00
MS52_scaled <- as.data.frame(scale(cfa_week0[1], scale = TRUE, center = TRUE))
MS52_raw <- cfa_week0[1]

MS00_raw <- cfa_week0[2]

MS01_raw <- cfa_week0[3]

source("https://bioconductor.org/biocLite.R")
biocLite("Rgraphviz")


#Scaling data of week 1

cfa1_week1 <- scale(cfa_week1, center = TRUE, scale = TRUE)
cfa1_week1 <- data.frame(cfa1_week1)
cfa1_week1 <- cbind(MS52_raw, cfa1_week1)
cfa1_week1

Others <- subset(Other_Variables, select=MAP:INFDOSE1)
Others <- as.data.frame(Others)
Categorical <- cbind(cfa1_week1, Other_Variables)


#Keep original factors of week 0 to see if factor 1 changes across time
model1 <- 'Factor 1 =~ HGB01 + HCT01 + RBC01 
           Factor 2 =~ STP01 + ALB01 + CAB01 + CHO01
Factor 3 =~ CK001 + ALT01 + AST01
Factor 4 =~ BUN01 + BC901 + BUA01 + P0401
Factor 5 =~ MCV01 + MCH01
Factor 6 =~ DDN01 + DBN01 + DJN01
AST01 ~~ 0*AST01
MCV01 ~~ 0*MCV01
DBN01 ~~ 0*DBN01
TLAMS52 ~ Factor 1 + Factor 2 + Factor 3 + Factor 4 + Factor 5 + Factor 6'

library(lavaan)
f1_week1 <- sem(model1, data = cfa1_week1, std.lv=TRUE, na.action=na.omit, missing = "ML")
summary(f1_week1, standardized = TRUE, fit.measures = TRUE)

library(semPlot)
semPaths(f1_week1,"std", whatLabels="std", intercepts=FALSE, style="lisrel", residuals = FALSE,
         nCharNodes=0, 
         nCharEdges=0,
         curveAdjacent = TRUE,title=TRUE, layout="tree2",curvePivot=FALSE)

p_week1 <- lavPredict(f1_week1)
p_week1 <- as.data.frame(p_week1)

colnames(p_week1) <- c("Factor1_w1", "Factor2_w1", "Factor3_w1", "Factor4_w1", "Factor5_w1", "Factor6_w1")

cfa2_week1 <- cbind(MS52_raw, p_week1)

factor_lm <- function(x) {
  f1 <- summary(lm(TLAMS52 ~ Factor1_w1, data = x))
  f2 <- summary(lm(TLAMS52 ~ Factor2_w1, data = x))
  f3 <- summary(lm(TLAMS52 ~ Factor3_w1, data = x))
  f4 <- summary(lm(TLAMS52 ~ Factor4_w1, data = x))
  f5 <- summary(lm(TLAMS52 ~ Factor5_w1, data = x))
  f6 <- summary(lm(TLAMS52 ~ Factor6_w1, data = x))
  return(list(f1,f2,f3,f4,f5,f6))
}

factor_lm(cfa2_week1)

#Scaling data for week 2
cfa1_week2 <- scale(cfa_week2, center = TRUE, scale = TRUE)
cfa1_week2 <- data.frame(cfa1_week2)
cfa1_week2 <- cbind(cfa1_week2, MS52_raw)
cfa1_week2

fa <- fa(cfa_week2, nfactors=6, rotate="varimax", fm="ml")
fa.diagram(fa)

#Keep original factors of week 0 to see if factor 1 changes across time

model2 <- 'Factor 1 =~ HGB02 + HCT02 + RBC02 
           Factor 2 =~ STP02 + ALB02 + CAB02 + CHO02
           Factor 3 =~ AST02 + CK002 + ALT02 
           Factor 4 =~ BUN02 + BC902 + BUA02 + P0402
           Factor 5 =~ MCV02 + MCH02
           Factor 6 =~ DBN02 + DDN02 + DJN02
           AST02 ~~ 0*AST02
           MCV02 ~~ 0*MCV02
           DBN02 ~~ 0*DBN02
           TLAMS52 ~ Factor 1 + Factor 2 + Factor 3 + Factor 4 + Factor 5 + Factor 6'

f1_week2 <- sem(model2, data = cfa1_week2, std.lv=TRUE, na.action=na.omit, missing = "fiml")
summary(f1_week2, standardized = TRUE, fit.measures = TRUE)


semPaths(f1_week2,"std", whatLabels="std", intercepts=FALSE, style="lisrel", residuals = FALSE,
         nCharNodes=0, 
         nCharEdges=0,
         curveAdjacent = TRUE,title=TRUE, layout="tree2",curvePivot=FALSE)

p_week2 <- lavPredict(f1_week2)
p_week2 <- as.data.frame(p_week2)

colnames(p_week2) <- c("Factor1_w2", "Factor2_w2", "Factor3_w2", "Factor4_w2", "Factor5_w2", "Factor6_w2")

cfa2_week2 <- cbind(MS52_raw, p_week2)

factor_lm <- function(x) {
  f1 <- summary(lm(TLAMS52 ~ Factor1_w2, data = x))
  f2 <- summary(lm(TLAMS52 ~ Factor2_w2, data = x))
  f3 <- summary(lm(TLAMS52 ~ Factor3_w2, data = x))
  f4 <- summary(lm(TLAMS52 ~ Factor4_w2, data = x))
  f5 <- summary(lm(TLAMS52 ~ Factor5_w2, data = x))
  f6 <- summary(lm(TLAMS52 ~ Factor6_w2, data = x))
  return(list(f1,f2,f3,f4,f5,f6))
}

factor_lm(cfa2_week2)

#scaling data for week 4
cfa1_week4 <- scale(cfa_week4, center = TRUE, scale = TRUE)
cfa1_week4 <- data.frame(cfa1_week4)
cfa1_week4


fa <- fa(cfa_week4[-1], nfactors=6, rotate="varimax", fm="ml")
fa.diagram(fa)

model4 <- 'Factor 1 =~ HGB04 + HCT04 + RBC04 
           Factor 2 =~ STP04 + ALB04 + CAB04 + CHO04
           Factor 3 =~ AST04 + CK004 + ALT04 
           Factor 4 =~ BUN04 + BC904 + BUA04 + P0404
           Factor 5 =~ MCV04 + MCH04
           Factor 6 =~ DBN04 + DDN04 + DJN04
           AST04 ~~ 0*AST04 
           BC904 ~~ 0*BC904
           MCH04 ~~ 0*MCH04 
           DBN04 ~~ 0*DBN04
           TLAMS52 ~ Factor 1 + Factor 2 + Factor 3 + Factor 4 + Factor 5 + Factor 6'

f1_week4 <- sem(model4, data = cfa1_week4, std.lv=TRUE, na.action=na.omit, missing = "fiml")
summary(f1_week4, standardized = TRUE, fit.measures = TRUE)

semPaths(f1_week4,"std", whatLabels="std", intercepts=FALSE, style="lisrel", residuals = FALSE,
         nCharNodes=0, 
         nCharEdges=0,
         curveAdjacent = TRUE,title=TRUE, layout="tree2",curvePivot=FALSE)

p_week4 <- lavPredict(f1_week4)
p_week4 <- as.data.frame(p_week4)

colnames(p_week4) <- c("Factor1_w4", "Factor2_w4", "Factor3_w4", "Factor4_w4", "Factor5_w4", "Factor6_w4")

cfa2_week4 <- cbind(cfa1_week4, p_week4)

factor_lm <- function(x) {
  f1 <- summary(lm(TLAMS52 ~ Factor1_w4, data = x))
  f2 <- summary(lm(TLAMS52 ~ Factor2_w4, data = x))
  f3 <- summary(lm(TLAMS52 ~ Factor3_w4, data = x))
  f4 <- summary(lm(TLAMS52 ~ Factor4_w4, data = x))
  f5 <- summary(lm(TLAMS52 ~ Factor5_w4, data = x))
  f6 <- summary(lm(TLAMS52 ~ Factor6_w4, data = x))
  return(list(f1,f2,f3,f4,f5,f6))
}

factor_lm(cfa2_week4)


#Factors for week 8
cfa1_week8 <- scale(cfa_week8, center = TRUE, scale = TRUE)
cfa1_week8 <- data.frame(cfa1_week8)
cfa1_week8 <- cbind(cfa1_week8, MS52_raw)
cfa1_week8


cor(cfa_week8, use = "pairwise.complete.obs")

fa <- fa(cfa_week8, nfactors=6, rotate="varimax", fm="ml")
fa.diagram(fa)

model8 <- 'Factor 1 =~ HGB08 + HCT08 + RBC08 
Factor 2 =~ STP08 + ALB08 + CAB08 + CHO08
Factor 3 =~ CK008 + ALT08 + AST08
Factor 4 =~ BUN08 + BC908 + BUA08 + P0408
Factor 5 =~ MCV08 + MCH08
Factor 6 =~ DDN08 + DBN08 + DJN08
AST08 ~~ 0*AST08
DBN08 ~~ 0*DBN08
BC908 ~~ 0*BC908
MCH08 ~~ 0*MCH08
TLAMS52 ~ Factor 1 + Factor 2 + Factor 3 + Factor 4 + Factor 5 + Factor 6'

library(lavaan)
f1_week8 <- sem(model8, data = cfa1_week8, std.lv=TRUE, na.action=na.omit, missing = "ML")
summary(f1_week8, standardized = TRUE, fit.measures = TRUE)

library(semPlot)
semPaths(f1_week8,"std", whatLabels="std", intercepts=FALSE, style="lisrel", residuals = FALSE,
         nCharNodes=0, 
         nCharEdges=0,
         curveAdjacent = TRUE,title=TRUE, layout="tree2",curvePivot=FALSE)

p_week8 <- lavPredict(f1_week8)
p_week8 <- as.data.frame(p_week8)

colnames(p_week8) <- c("Factor1_w8", "Factor2_w8", "Factor3_w8", "Factor4_w8", "Factor5_w8", "Factor6_w8")

#Factors for week 52
cfa1_week52 <- scale(cfa_week52, center = TRUE, scale = TRUE)
cfa1_week52 <- data.frame(cfa1_week52)
cfa1_week52 <- cbind(cfa1_week52, MS52_raw)
cfa1_week52

fa <- fa(cfa_week52, nfactors=6, rotate="varimax", fm="ml")
fa.diagram(fa)

model52 <- 'Factor 1 =~ HGB52 + HCT52 + RBC52 
Factor 2 =~ STP52 + ALB52 + CAB52 + CHO52
Factor 3 =~ CK052 + ALT52 + AST52
Factor 4 =~ BUN52 + BC952 + BUA52 + P0452
Factor 5 =~ MCV52 + MCH52
Factor 6 =~ DDN52 + DBN52 + DJN52
AST52 ~~ 0*AST52
DBN52 ~~ 0*DBN52
BC952 ~~ 0*BC952
HCT52 ~~ 0*HCT52
MCV52 ~~ 0*MCV52
TLAMS52 ~ Factor 1 + Factor 2 + Factor 3 + Factor 4 + Factor 5 + Factor 6'

library(lavaan)
f1_week52 <- sem(model52, data = cfa1_week52, std.lv=TRUE, na.action=na.omit, missing = "ML")
summary(f1_week52, standardized = TRUE, fit.measures = TRUE)

library(semPlot)
semPaths(f1_week52,"std", whatLabels="std", intercepts=FALSE, style="lisrel", residuals = FALSE,
         nCharNodes=0, 
         nCharEdges=0,
         curveAdjacent = TRUE,title=TRUE, layout="tree2",curvePivot=FALSE)

p_week52 <- lavPredict(f1_week52)
p_week52 <- as.data.frame(p_week52)

colnames(p_week52) <- c("Factor1_w52", "Factor2_w52", "Factor3_w52", "Factor4_w52", "Factor5_w52", "Factor6_w52")

#Combined dataframe
Factors_Overtime <- do.call("cbind", list(MS52_raw, p_week0, p_week1, p_week2, p_week4, p_week8, p_week52, Other_Variables))

Factors_Overtime$SEXCD <- as.factor(Factors_Overtime$SEXCD) 
Factors_Overtime$AGE <- as.numeric(Factors_Overtime$AGE)

#UCP is not working for long format so I used the wide format for each timepoint
#Week 0
Factor1 <- cfa_week0[, c("RBC00", "HGB00", "HCT00")]
Factors_W0 <- do.call("cbind", list(MS52_raw, p_week0, Other_Variables, Factor1))

#Turning variables into factor
Factors_W0 <- Factors_W0 %>% mutate_if(is.character,as.factor) #They only take factor values
Factors_W0$SEXCD <- as.factor(Factors_W0$SEXCD) 
Factors_W0$DEATHRPD <- as.factor(Factors_W0$DEATHRPD)
Factors_W0$AGE <- as.numeric(Factors_W0$AGE)

#ctree can't take missing values in response!
Factors_W0 <- subset(Factors_W0, !is.na(TLAMS52))

#UCP for LMS52
library(party)
UCP_Factors_W0 <- ctree(TLAMS52 ~ Factor1_w0 + AGE + SPLVL1 + ASIMPC01_A, data = Factors_W0)
plot(UCP_Factors_W0)


#UCP for week 4
Factors_W4 <- do.call("cbind", list(MS52_raw, p_week4, Other_Variables, Factor1))

Factors_W4 <- Factors_W4 %>% mutate_if(is.character,as.factor) #They only take factor values
Factors_W4$SEXCD <- as.factor(Factors_W4$SEXCD) 

Factors_W4 <- subset(Factors_W4, !is.na(TLAMS52))

UCP_Factors_W4 <- ctree(TLAMS52 ~ Factor1 + Factor2 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_W4)
plot(UCP_Factors_W4)


#UCP for week1-4
Factors_Overtime <- Factors_Overtime %>% mutate_if(is.character,as.factor)
Factors_Overtime$SEXCD <- as.factor(Factors_Overtime$SEXCD) 

Factors_Overtime <- subset(Factors_Overtime, !is.na(TLAMS52))

UCP_Factor1_Overtime <- ctree(TLAMS52 ~ Factor1_w0 + Factor1_w1 + Factor1_w2 + Factor1_w4 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor1_Overtime)

UCP_Factor2_Overtime <- ctree(TLAMS52 ~ Factor2_w0 + Factor2_w1 + Factor2_w2 + Factor2_w4 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor2_Overtime)

UCP_Factor3_Overtime <- ctree(TLAMS52 ~ Factor3_w0 + Factor3_w1 + Factor3_w2 + Factor3_w4 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor3_Overtime)

UCP_Factor4_Overtime <- ctree(TLAMS52 ~ Factor4_w0 + Factor4_w1 + Factor4_w2 + Factor4_w4 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor4_Overtime)

UCP_Factor5_Overtime <- ctree(TLAMS52 ~ Factor5_w0 + Factor5_w1 + Factor5_w2 + Factor5_w4 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor5_Overtime)

UCP_Factor6_Overtime <- ctree(TLAMS52 ~ Factor6_w0 + Factor6_w1 + Factor6_w2 + Factor6_w4 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor6_Overtime)

#Transforming data to long format

#Taking Correlation r values between MS52 and Factors across time and put into df
Cor_Factors <- cor(Factors_Overtime[1:37], use="pairwise.complete.obs")
Cor_Factors <- head(Cor_Factors, 1)
Cor_Factors <- as.data.frame(Cor_Factors)

#Transform into long format
library(reshape)
Cor_Factors_long <- melt(Cor_Factors, id.vars=c("TLAMS52"),
                         measure.vars=c("Factor1_w0", "Factor1_w1", "Factor1_w2", "Factor1_w4", "Factor1_w8", "Factor1_w52",
                                        "Factor2_w0", "Factor2_w1", "Factor2_w2", "Factor2_w4", "Factor2_w8", "Factor2_w52",
                                        "Factor3_w0", "Factor3_w1", "Factor3_w2", "Factor3_w4", "Factor3_w8", "Factor3_w52",
                                        "Factor4_w0", "Factor4_w1", "Factor4_w2", "Factor4_w4", "Factor4_w8", "Factor4_w52",
                                        "Factor5_w0", "Factor5_w1", "Factor5_w2", "Factor5_w4", "Factor5_w8", "Factor5_w52",
                                        "Factor6_w0", "Factor6_w1", "Factor6_w2", "Factor6_w4", "Factor6_w8", "Factor6_w52"),
                         variable.name ="Time After Injury", 
                         value.name ="Factors")

#Splitting into 2 separate columns
library(tidyr)
Cor_Factors_long <- separate(data = Cor_Factors_long, col = variable, 
                             into = c("Factors", "TimeAfterInjury"), sep = "\\_")

Cor_Factors_long$TimeAfterInjury = factor(Cor_Factors_long$TimeAfterInjury, levels=c('w0','w1','w2','w4','w8','w52'))

#Plot across time
library(ggplot2)
p <- ggplot(Cor_Factors_long, aes(x=TimeAfterInjury, y=value)) +
  geom_jitter()+
  geom_smooth(color='red', method='lm', formula=y~x, aes(group=1), se=FALSE) +
  ggtitle("Correlation between Motor Scores at week 52 and Factors from week 0 to week 52")+
  theme_bw() +
  labs(x="Time After Injury", y="r values at MS week 52") +
  facet_grid(.~Factors)
p

#UCP for week0-52
library(dplyr)
Factors_Overtime <- Factors_Overtime %>% mutate_if(is.character,as.factor)



Factors_Overtime <- subset(Factors_Overtime, !is.na(TLAMS52))

UCP_Factor1_Overtime <- ctree(TLAMS52 ~ Factor1_w0 + Factor1_w1 + Factor1_w2 + Factor1_w4 + Factor1_w8 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor1_Overtime)

UCP_Factor2_Overtime <- ctree(TLAMS52 ~ Factor2_w0 + Factor2_w1 + Factor2_w2 + Factor2_w4 + Factor2_w8 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor2_Overtime)

UCP_Factor3_Overtime <- ctree(TLAMS52 ~ Factor3_w0 + Factor3_w1 + Factor3_w2 + Factor3_w4 + Factor3_w8 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor3_Overtime)

UCP_Factor4_Overtime <- ctree(TLAMS52 ~ Factor4_w0 + Factor4_w1 + Factor4_w2 + Factor4_w4 + Factor4_w8 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor4_Overtime)

UCP_Factor5_Overtime <- ctree(TLAMS52 ~ Factor5_w0 + Factor5_w1 + Factor5_w2 + Factor5_w4 + Factor5_w8 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor5_Overtime)

UCP_Factor6_Overtime <- ctree(TLAMS52 ~ Factor6_w0 + Factor6_w1 + Factor6_w2 + Factor6_w4 + Factor6_w8 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factor6_Overtime)

#Tree for individual week!
UCP_week0 <- ctree(TLAMS52 ~ Factor1_w0 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_week0)

UCP_week1 <- ctree(TLAMS52 ~ Factor1_w1 + Factor2_w1 +  Factor6_w1 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = subset(Factors_Overtime, !is.na(TLAMS52)))
plot(UCP_week1)

UCP_week2 <- ctree(TLAMS52 ~ Factor1_w2 + Factor2_w2 +Factor3_w2 + Factor4_w2 + Factor5_w2 + Factor6_w2 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_week2)

UCP_week4 <- ctree(TLAMS52 ~ Factor1_w4 + Factor2_w4 +Factor3_w4 + Factor4_w4 + Factor5_w4 + Factor6_w4 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_week4)

UCP_week8 <- ctree(TLAMS52 ~ Factor1_w8 + Factor2_w8 +Factor3_w8 + Factor4_w8 + Factor5_w8 + Factor6_w8 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_week8)

UCP_week52 <- ctree(TLAMS52 ~ Factor1_w52 + Factor2_w52 +Factor3_w52 + Factor4_w52 + Factor5_w52 + Factor6_w52 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_week52)

#Models for week 0! (Factor 1 is strong)
library(mixlm)
library(MASS)
drop1(lm_week0_factor1 , test = "Chi")
step <- stepAIC(lm_week0_all, direction = "both")
summary(step)
forward(lm_week0_all)
backward(lm_week0_all)

lm_week0_all <- lm(TLAMS52 ~ Factor1_w0 + Factor2_w0 + Factor3_w0 + Factor4_w0 + Factor5_w0 + Factor6_w0 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)

lm_week0_factor1 <- lm(TLAMS52 ~  Factor1_w0, data = Factors_Overtime)
lm_week0_factor2 <- lm(TLAMS52 ~  Factor2_w0, data = Factors_Overtime)
lm_week0_factor3 <- lm(TLAMS52 ~  Factor3_w0, data = Factors_Overtime)
lm_week0_factor4 <- lm(TLAMS52 ~  Factor4_w0, data = Factors_Overtime)
lm_week0_factor5 <- lm(TLAMS52 ~  Factor5_w0, data = Factors_Overtime)
lm_week0_factor6 <- lm(TLAMS52 ~  Factor6_w0, data = Factors_Overtime)

lm_week0_factors <- lm(TLAMS52 ~  Factor1_w0 + Factor2_w0 + Factor3_w0, data = Factors_Overtime)
lm_week0_inj <- lm(TLAMS52 ~ SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)


summary(lm_week0_factor1)
summary(lm_week0_factor2)
summary(lm_week0_factor3)
summary(lm_week0_factors)
summary(lm_week0_inj)
summary(lm_week0_all)

AIC(lm_week0_factor1)
AIC(lm_week0_factors)
AIC(lm_week0_inj)

Factors_Overtime <- subset(Factors_Overtime, !is.na(Factor1_w0))
Factors_Overtime <- subset(Factors_Overtime, !is.na(ASIMPC01_A))

anova(lm_week0_factor1, lm_week0_factors, lm_week0_inj)
anova(lm_week0_factor1, lm_week0_inj)

qqnorm(resid(lm_week0))
qqline(resid(lm_week0))

plot(fitted(lm_week0_inj), resid(lm_week0_inj))
abline(a=0,b =0)

#Model diagnostic
library(car)
crPlot(lm_week0_factor1, "Factor1_w0")
vif(lm_week0_factor1)

#Models for week 1 (Factor 1 & 2 is strong) 
drop1(lm_week1_factor1, test = "Chi")
lm_week1_all <- lm(TLAMS52 ~ Factor1_w1 + Factor2_w1 + Factor3_w1 + Factor4_w1 + Factor5_w1 + Factor6_w1 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)

step <- stepAIC(lm_week1_factors, direction = "both")
summary(step)
forward(lm_week1_factors)
backward(lm_week1_factors)

lm_week1_factor1and2 <- lm(TLAMS52 ~ Factor1_w1 + Factor2_w1 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
lm_week1_factor1 <- lm(TLAMS52 ~  Factor1_w1 + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)
lm_week1_factor2 <- lm(TLAMS52 ~  Factor2_w1 + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)
lm_week1_factor3 <- lm(TLAMS52 ~  Factor3_w1 + SEXCD + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)

autoplot(lm_week1_factor2, which = 1:6, colour = 'ASIMPC01_A')


lm_week1_factors <- lm(TLAMS52 ~ Factor1_w1 + Factor2_w1 + Factor3_w1 + Factor4_w1 + Factor5_w1 + Factor6_w1, data = Factors_Overtime)
lm_week1_inj <- lm(TLAMS52 ~ SEXCD + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)

summary(lm_week1_factor1and2)
summary(lm_week1_factor1)
summary(lm_week1_factor2)
summary(lm_week1_factor3)
summary(lm_week1_factors)
summary(lm_week1_inj)

Factors_Overtime <- subset(Factors_Overtime, !is.na(ASIMPC01_A))
Factors_Overtime <- subset(Factors_Overtime, !is.na(Factor1_w1))



AIC(lm_week1_factor1)
AIC(lm_week1_factor2)
AIC(lm_week1_factors)
AIC(lm_week1_inj)
anova(lm_week1_factor1, lm_week1_factors, lm_week1_inj, test = "Chisq")
anova(lm_week1_factor1, lm_week1_inj, test = "Chisq")
anova(lm_week1_factor1, lm_week1_factor1and2, test = "Chisq")

#Model Diagnostics 
crPlots(lm_week1_factor2)
vif(lm_week1_factor1)

#Models for week 2 (Factor 1 & 2, 6 is strong)
lm_week2_all <- lm(TLAMS52 ~ Factor1_w2 + Factor2_w2 + Factor3_w2 + Factor4_w2 + Factor5_w2 + Factor6_w2, data = Factors_Overtime)

step <- stepAIC(lm_week2_all, direction = "both")
summary(step)
forward(lm_week2_all)
backward(lm_week2_all)

lm_week2_factor2and6 <- lm(TLAMS52 ~ Factor2_w2 + Factor6_w2 + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)

lm_week2_factor1 <- lm(TLAMS52 ~  Factor1_w2 + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)
lm_week2_factor2 <- lm(TLAMS52 ~  Factor2_w2 + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)
lm_week2_factor6 <- lm(TLAMS52 ~  Factor6_w2 + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)
lm_week2_factor3 <- lm(TLAMS52 ~  Factor3_w2 + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)

lm_week2_factors <- lm(TLAMS52 ~  Factor1_w2 + Factor2_w2 + Factor6_w2, data = Factors_Overtime)
lm_week2_inj <- lm(TLAMS52 ~ AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)

summary(lm_week2_factor2and6)
summary(lm_week2_factor1)
summary(lm_week2_factor2)
summary(lm_week2_factor6)
summary(lm_week2_factor3)
summary(lm_week2_factors)
summary(lm_week2_inj)

Factors_Overtime <- subset(Factors_Overtime, !is.na(ASIMPC01_A))
Factors_Overtime <- subset(Factors_Overtime, !is.na(Factor1_w2))
Factors_Overtime <- subset(Factors_Overtime, !is.na(Factor2_w2))


AIC(lm_week2_factor1)
AIC(lm_week2_factor2)
AIC(lm_week2_factors)
AIC(lm_week2_inj)

anova(lm_week2_factor6, lm_week2_factors, test = "Chisq")
anova(lm_week2_factor6, lm_week2_inj, test = "Chisq")
anova(lm_week2_factor1, lm_week2_factor2, lm_week2_factor6)

#model diagnostic 
crPlot(lm_week2_factor2, "Factor1_w2")
vif(lm_week2_factor1)

#Models for week 4 (factor1, 2, and 6)
drop1(lm_week4, test = "Chi")

step <- stepAIC(lm_week4_factors, direction = "both")
summary(step)
forward(lm_week4_inj)
backward(lm_week4_inj)

lm_week4_factor1 <- lm(TLAMS52 ~  Factor1_w4 + SEXCD + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)
lm_week4_factor2 <- lm(TLAMS52 ~  Factor2_w4 + SEXCD + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)
lm_week4_factor3 <- lm(TLAMS52 ~  Factor3_w4 + SEXCD + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)
lm_week4_factor4 <- lm(TLAMS52 ~  Factor4_w4 + SEXCD + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)
lm_week4_factor5 <- lm(TLAMS52 ~  Factor5_w4 + SEXCD + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)
lm_week4_factor6 <- lm(TLAMS52 ~  Factor6_w4 + SEXCD + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)


lm_week4_factors <- lm(TLAMS52 ~  Factor1_w4 + Factor2_w4 + Factor3_w4 + Factor4_w4 + Factor5_w4 + Factor6_w4, data = Factors_Overtime)
lm_week4_inj <- lm(TLAMS52 ~ SEXCD + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime)

summary(lm_week4_factor1)
summary(lm_week4_factor2)
summary(lm_week4_factor3)
summary(lm_week4_factor4)
summary(lm_week4_factor5)
summary(lm_week4_factor6)

summary(lm_week4_factors)
summary(lm_week4_inj)

Factors_Overtime <- subset(Factors_Overtime, !is.na(ASIMPC01_A))
Factors_Overtime <- subset(Factors_Overtime, !is.na(Factor1_w2))



AIC(lm_week4_factor1)
AIC(lm_week4_factor2)
AIC(lm_week4_factor6)
AIC(lm_week4_factors)
AIC(lm_week4_inj)
anova(lm_week4_factor1, lm_week4_factors, lm_week4_inj, test = "Chisq")
anova(lm_week4_factor1, lm_week4_inj, test = "Chisq")
anova(lm_week4_factor1, lm_week4_factor2, lm_week4_factor6, test = "Chisq")
anova(lm_week4_factor1, lm_week4_factor6)

#model diagnostic
crPlots(lm_week4_factor6)
vif(lm_week4_factor1)

#Identifying Anemia 
cfa_week0$Anemia <- ifelse(cfa_week0$HGB00 <= 12.3, 1, 0)
cfa_week0$Anemia <- as.factor(cfa_week0$Anemia)
Anemia <- as.data.frame(cfa_week0$Anemia)
colnames(Anemia) <- c("Anemia")

ggplot(data = na.omit(cfa_week0), aes(x=Anemia)) +
  geom_bar(width = 0.3) +
  theme_bw() + 
  ggtitle("Numbers of patients having Anemia at week 0")

cfa_week0$Type_Anemia <- ifelse(cfa_week0$Anemia == 0, "0",
                                ifelse(cfa_week0$MCV00 < 80, "1",
                                       ifelse(cfa_week0$MCV00 > 100, "3", "2")))


ggplot(data = na.omit(cfa_week0), aes(x=Type_Anemia)) +
  geom_bar(width = 0.3) +
  theme_bw() + 
  ggtitle("Numbers of patients having different types of Anemia at week 0")

#Is Anemia due to excessive blood loss?
Other_Variables$Pulse <- Other_Variables$BPSYS - Other_Variables$BPDIA

#Is Anemia having an impact on MS52?

#Combine dataframe in week 0 to identify if Anemia matters

Factors_Overtime_Anemia <- cbind(Factors_Overtime, Anemia)

Factors_Overtime_Anemia <- subset(Factors_Overtime_Anemia, !is.na(ASIMPC01_A))

Factors_Overtime_Anemia <- subset(Factors_Overtime_Anemia, !is.na(Anemia))
Factors_Overtime_Anemia <- subset(Factors_Overtime_Anemia, !is.na(TLAMS52))

Factors_Overtime_Anemia <- Factors_Overtime_Anemia %>% mutate_if(is.character,as.factor)

Factors_Overtime_Anemia$SEXCD <- as.factor(Factors_Overtime_Anemia$SEXCD) 
Factors_Overtime_Anemia$AGE <- as.numeric(Factors_Overtime_Anemia$AGE)
Factors_Overtime_Anemia$DEATHRPD <- as.factor(Factors_Overtime_Anemia$DEATHRPD)
Factors_Overtime_Anemia$Factor1_w0_flipped <- Factors_Overtime_Anemia$Factor1_w0+1

lm_week0_Anemia <- lm(TLAMS52 ~ Factor1_w0_flipped +  AGE + SEXCD  + ASIMPC01_A + SPLVL1, data = Factors_Overtime_Anemia)
summary(lm_week0_Anemia)
sjp.lm(lm_week0_Anemia)
drop1(lm_week0_Anemia, test = "Chi")

lm_week0_Anemia_inj <- lm(TLAMS52 ~ SEXCD + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Anemia)
summary(lm_week0_Anemia_inj)

lm_week0_Anemia_alone <- lm(TLAMS52 ~ Anemia, data = Factors_Overtime_Anemia)
summary(lm_week0_Anemia_alone)

lm_week0_Anemia_Factor <- lm(TLAMS52 ~ Factor1_w0_flipped, data = Factors_Overtime_Anemia)
summary(lm_week0_Anemia_Factor)


AIC(lm_week0_Anemia)
AIC(lm_week0_Anemia_inj)
AIC(lm_week0_Anemia_alone)
AIC(lm_week0_Anemia_Factor)

anova(lm_week0_Anemia, lm_week0_Anemia_Factor, test = "Chisq")

#Diagnostic
qqnorm(resid(lm_week0_Anemia))
qqline(resid(lm_week0_Anemia))

#Diagnostic plot -> looks like we have heteroscedascity and non-normality problem
autoplot(lm_week0_Anemia, which = 1:6, colour = 'ASIMPC01_A')
plot(lm_week0_Anemia, add.smooth = TRUE)


#Skewed data
library(moments)
skewness(lm_week0_Anemia_Factor$residuals)

plot(Effect(c("Factor1_w0", "ASIMPC01_A"), lm_week0_Anemia_Factor), multiline = TRUE)

#Logarithm is clumsy and doesn't work
Factors_Overtime_Anemia$log.MS52 <- log(Factors_Overtime_Anemia$TLAMS52 + 1)
Factors_Overtime_Anemia$log.Factor1_w0 <- log(Factors_Overtime_Anemia$Factor1_w0 + 1)

#Graph of lm and diagnostic 
library(sjPlot)
sjp.lm(lm_week0_Anemia)

#Model Diagnostic as follows: Outliers -> Linearity -> Heteroscedasticiy -> Non-normality

#If we know ahead of time of the dist (skewed), can use non-parametric test 
#Parametric follows certain dist (uses means)), non-paramentric doesn't follow (uses medians)

#Generalized Additive Models <- semipametric!
# Can be used to address non-linearity and (heteroscedasticity?)
library(mgcv)
gam_week0_Anemia_Factor <- gam(TLAMS52 ~ Factor1_w0_flipped + AGE + SPLVL1 + ASIMPC01_A + SEXCD, data = Factors_Overtime_Anemia)
summary(gam_week0_Anemia_Factor)


#Outliers
outlierTest(lm_week0_Anemia)
cooksd <- cooks.distance(lm_week0_Anemia)
plot(cooksd, pch="*", cex=1, main="Influential Obs by Cooks distance")
abline(h = 4*mean(cooksd, na.rm=T), col="red") 
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")

#Best Outlier test!!! I think - no outliers
source("http://goo.gl/UUyEzD")
b <- outlierKD(Factors_Overtime, Factor1_w0)

b <- cbind(Factors_Overtime_Anemia, b)
#No outliers worth removing 

#Or using robust linear regression to addressing for outliers! (and heteroscedascity?)
rlm_week0_Anemia_Factor <- rlm(TLAMS52~ Factor1_w0_flipped + AGE + SPLVL1 + ASIMPC01_A + SEXCD, data = Factors_Overtime_Anemia, init = "ls", maxit = 50)

summary(rlm_week0_Anemia_Factor)
library(sfsmisc)
f.robftest(rlm_week0_Anemia_Factor, var = "Factor1_w0_flipped")

#get the outlier observation
cDist <- cooks.distance (lm_week0_Anemia_Factor) # get cooks distance

n <- nrow (Factors_Overtime_Anemia) # No of rows

resids <- stdres (lm_week0_Anemia_Factor) #errors

a <- cbind(Factors_Overtime_Anemia, cDist, resids)

cDist > 4/797

#Check for linearity - Yes
crPlots(lm_week0_Anemia)

#Multicolinearity - No (Cut off is 2, 1 is no colinearity)
library(car)
vif(lm_week0_Anemia)

#Other tests
library(gvlma)
gvmodel <- gvlma(lm_week0_Anemia) 
summary(gvmodel)

#Normal distribution? -> not
shapiro.test(resid(lm_week0_Anemia))

#Heteroscasdicity -> yes!
library(lmtest)
ncvTest(lm_week0_Anemia)
bptest(lm_week0_Anemia)

spreadLevelPlot(lm_week0_Anemia)
# We can also try BoxCox but with our dataset, boxcox may not be appropriate because of extreme values (0)

#Weighted regression <- don't understand!!!!!
wl = 1 / apply(Factors_Overtime_Anemia, 2, function(x){ var(x, na.rm=T) })
Factors_Overtime_Anemia$w = with(Factors_Overtime_Anemia, ifelse(ASIMPC01_A=="A", wl[3], 
                                           ifelse(ASIMPC01_A=="B", wl[3], wl[3])))

Factors_Overtime_Anemia$Factor1_w0_adjusted <- Factors_Overtime_Anemia$Factor1_w0 + 3

wts <- 1/fitted(lm(abs(residuals(lm_week0_Anemia_Factor)) ~ Factor1_w0))^2

lm_week0_Anemia_Factor_weighted <- lm(TLAMS52 ~ Factor1_w0_adjusted + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Anemia, weights = Factor1_w0_adjusted)
summary(lm_week0_Anemia_Factor_weighted)

plot(lm_week0_Anemia_Factor_weighted)
bptest(lm_week0_Anemia_Factor_weighted)
#I believe weights should be on ASIMPC01_A since they are the one driving heteroscasdicity

# Option A for BoxCox transformation
library(caret)
boxcoxMS <- BoxCoxTrans(Factors_Overtime_Anemia$TLAMS52+1)
print(boxcoxMS)
new_boxcox <- cbind(Factors_Overtime_Anemia, MS52_boxcox=predict(boxcoxMS, Factors_Overtime_Anemia$TLAMS52+1))
lm_new_boxcox <- lm(MS52_boxcox ~  Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = new_boxcox)
summary(lm_new_boxcox)
ncvTest(lm_new_boxcox)
bptest(lm_new_boxcox)
shapiro.test(resid(lm_new_boxcox))
qqnorm(resid(lm_new_boxcox))
qqline(resid(lm_new_boxcox))
hist(new_boxcox$MS52_boxcox)

#Option B for BoxCox transformation
boxcox_lamda <- boxcox(TLAMS52+1 ~  Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
lambda <- boxcox_lamda$x[which.max(boxcox_lamda$y)]
lm_boxcox <- lm((((TLAMS52+1)^lambda - 1)/lambda) ~ Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
summary(lm_boxcox)
ncvTest(lm_boxcox)
bptest(lm_boxcox)
shapiro.test(resid(lm_boxcox))
qqnorm(resid(lm_boxcox))

#Option C for BoxCox transformation 
summary(p1 <- powerTransform(lm_week0_Anemia))

# => Both options did NOT transform data!!!!!!!!
# ncvTest and bptest are different because bptest is studendized and ncvTest is not
# => which means bptest is more robust than the original one (ncvTest)
# => go with bptest -> null is: homoscedasicity 

#Use this to treat for heteroscadasicity or generalized additive measures or rlm?
library(lmtest)
library(sandwich)
library(caret)
vcovHC(lm_week0_Anemia)
coeftest(lm_week0_Anemia, vcov. = vcovHC)

coeftest(lm_week0_Anemia,vcov=hccm(lm_week0_Anemia))
coeftest(lm_week0_Anemia, vcov. = vcovHC(lm_week0_Anemia)) # same as above

HC <- hccm(lm_week0_Anemia) 
homo <- coeftest(lm_week0_Anemia, vcov. = HC) #Same as above too

waldtest(lm_week0_Anemia,  vcov = vcovHC)

#Look at histogram of response
hist(new_boxcox$MS52_boxcox)

#ctree for Anemia
UCP_Anemia <- ctree(TLAMS52 ~ Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
plot(UCP_Anemia)

str(Factors_Overtime_Anemia)

ex <- lm(TLAMS52 ~ ASIMPC01_A*SPLVL1 + ASIMPC01_A*Factor1_w0 +  Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
summary(ex)

#UCP for individual component of factor 1
RBC_columns <- c("RBC00", "HCT00", "HGB00")
RBC <- cfa_week0[RBC_columns]
RBC_scaled <- scale(RBC, scale = TRUE, center = TRUE)
RBC_scaled <- as.data.frame(RBC_scaled)

Factors_Overtime_Anemia <- cbind(Factors_Overtime_Anemia, RBC_scaled)

UCP_RBC <- partykit::ctree(TLAMS52 ~ RBC00 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
plot(UCP_RBC)

UCP_HGB <- partykit::ctree(TLAMS52 ~ HGB00 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
plot(UCP_HGB)

UCP_HCT <- partykit::ctree(TLAMS52 ~ HCT00 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
plot(UCP_HCT)

#Figures for Factor 1
ggplot(data= subset(Factors_Overtime_Anemia, !is.na(ASIMPC01_A)), aes(x=TLAMS52, y=Factor1_w0_flipped))+ 
  geom_jitter(shape=1) + 
  geom_smooth(method = lm) + 
  scale_colour_hue(l=50) +
  theme_bw() +
  ggtitle("Factor 1 (Anemia) Across ASI in Motor Scores") + 
  facet_grid(. ~ ASIMPC01_A, scales = "free_x")

#Walking functions 
Other_Variables$Walk <- ifelse(Other_Variables$MODBEN52 == 5, "1",
                               ifelse(Other_Variables$MODBEN52 == 6, "1",
                                      ifelse(Other_Variables$MODBEN52 == 7, "1", "0")))
Other_Variables$Walk <- as.factor(Other_Variables$Walk)

Other_Variables <- Other_Variables %>% mutate_if(is.character,as.factor)

Factors_Walk <- do.call("cbind", list(MS52_raw, p_week0, p_week1, p_week2, p_week4, p_week8, p_week52, Other_Variables))

Factors_Walk$SEXCD <- as.factor(Factors_Walk$SEXCD) 
Factors_Walk$AGE <- as.numeric(Factors_Walk$AGE)
Factors_Walk$ASIMPC01_A <- as.factor(Factors_Walk$ASIMPC01_A)



glm_walk <- glm(Walk ~ Factor1_w0 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Walk, family = "binomial")
summary(glm_walk)

drop1(glm_walk, test = "Chi")
stepAIC(glm_walk)

glm_walk_reduce <- glm(Walk ~ Factor1_w0 + AGE + ASIMPC01_A, data = Factors_Walk, family = "quasibinomial")
summary(glm_walk_reduce)

anova(glm_walk, glm_walk_reduce, test = "Chisq")

Factors_Walk <- subset(Factors_Walk, !is.na(Walk))

UCP_walk <- ctree(Walk ~ Factor1_w0 + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Walk)

plot(UCP_walk)

sjp.glm(glm_walk_reduce)
#Test for dispersion in glm model
library(msme)
P__disp(glm_walk)

#Marked Recovery 
# 9 in MODBEN52 means NA so I transformed them
MODBEN52 <- as.data.frame(Other_Variables$MODBEN52)
MODBEN52[MODBEN52 == 9] <- NA
colnames(MODBEN52) <- c("MODBEN52_1")

Other_Variables <- cbind(Other_Variables, MODBEN52)

#Turning ASI grades into numeric 
Other_Variables$ASI_numeric <- ifelse(Other_Variables$ASIMPC01_A == "A", 1,
                                          ifelse(Other_Variables$ASIMPC01_A == "B", 2,
                                                 ifelse(Other_Variables$ASIMPC01_A == "C", 3, 4)))

#Calculate Marked Recovery 
Other_Variables$Marked_Recovery <- ifelse(Other_Variables$MODBEN52_1 - Other_Variables$ASI_numeric >= 2, 1,0)


Factors_MR <- do.call("cbind", list(MS01_raw, MS00_raw, MS52_raw, p_week0, p_week1, p_week2, p_week4, p_week8, p_week52, Other_Variables, RBC))

# GLM for marked Recovery (Factor 1 and 2 alone are significant but when consider other variables, only Factor 1 is sig)
Factors_MR <- Factors_MR %>% mutate_if(is.character,as.factor)
Factors_MR$SEXCD <- as.factor(Factors_MR$SEXCD) 
Factors_MR$AGE <- as.numeric(Factors_MR$AGE)
Factors_MR$Marked_Recovery <- as.factor(Factors_MR$Marked_Recovery)
Factors_MR$Factor1_w0_flipped <- Factors_MR$Factor1_w0+1

glm_MR_factor1 <- glm(Marked_Recovery ~ Factor1_w0, data = Factors_MR, family = quasibinomial)
summary(glm_MR_factor1)

glm_MR_factor2 <- glm(Marked_Recovery ~ Factor2_w0, data = Factors_MR, family = quasibinomial)

glm_MR_factor3 <- glm(Marked_Recovery ~ Factor3_w0, data = Factors_MR, family = quasibinomial)

glm_MR_factor4 <- glm(Marked_Recovery ~ Factor4_w0, data = Factors_MR, family = quasibinomial)

glm_MR_factor5 <- glm(Marked_Recovery ~ Factor5_w0, data = Factors_MR, family = quasibinomial)

glm_MR_factor6 <- glm(Marked_Recovery ~ Factor6_w0, data = Factors_MR, family = quasibinomial)


glm_MR_inj <- glm(Marked_Recovery ~ SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_MR, family = quasibinomial)
summary(glm_MR_inj)

anova(glm_MR_factor1, glm_MR_reduced, test = "Chisq")

glm_MR <- glm(Marked_Recovery ~ Factor1_w0_flipped +  AGE + SEXCD  + ASIMPC01_A + SPLVL1, data = Factors_MR, family = binomial)
summary(glm_MR)

glm_MR_2 <- glm(Marked_Recovery ~ Factor2_w0 +  AGE + SEXCD  + ASIMPC01_A + SPLVL1, data = Factors_MR, family = quasibinomial)
summary(glm_MR_2)

drop1(glm_MR_inj, test = "Chi")

glm_MR_reduced <- glm(Marked_Recovery ~ Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_MR, family = quasibinomial)
summary(glm_MR_reduced)


anova(glm_MR_reduced, glm_MR, test = "Chisq")

#test for dispersion 
#One way to see if model is overdispersed is looking at Residual Deviance : degrees of freedom
# more than 1 -> Overdispersion
library(msme)
P__disp(glm_MR)

#Ctree
Factors_MR <- subset(Factors_MR, !is.na(Marked_Recovery))

UCP_MR <- ctree(Marked_Recovery ~ Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_MR)
plot(UCP_MR)

#Model diagnostics - Test for linearity and multicollinearity 
crPlots(glm_MR)
vif(glm_MR)

library(boot)
glm.diag.plots (glm_MR)


library(pROC) #Only works with numeric or ordered factors so if it is factors, have to order it
Factors_MR$Ordered_ASIMPC01_A <- factor(Factors_MR$ASIMPC01_A, levels=c("A", "B", "C", "D"), ordered=TRUE)
Factors_MR$Ordered_SPLVL1 <- factor(Factors_MR$SPLVL1, levels=c("T", "C"), ordered= TRUE)
Factors_MR$Ordered_SEXCD <- factor(Factors_MR$SEXCD, levels=c("1", "2"), ordered= TRUE)
Factors_MR$Ordered_Factor1_w0 <- factor(Factors_MR$Factor1_w0, ordered = TRUE)
AUC = roc(Marked_Recovery ~ Factor1_w0 + Ordered_ASIMPC01_A + Ordered_SPLVL1 + AGE + Ordered_SEXCD, data=Factors_MR, plot = TRUE, percent = TRUE, ci = TRUE) 

F1 = roc(Marked_Recovery ~ Ordered_Factor1_w0, data=Factors_MR) 
plot(F1, col="red", main = "ROC for Factor 1 (Anemia)")

F2 = roc(Marked_Recovery ~ Ordered_ASIMPC01_A, data=Factors_MR) 
plot(F2, col="red", main = "ROC for AIS grades")

F3 = roc(Marked_Recovery ~ Ordered_SPLVL1, data=Factors_MR) 
plot(F3, col="red", main = "ROC for Spinal Levels")

F4 = roc(Marked_Recovery ~ AGE, data=Factors_MR) 
plot(F4, col="red", main = "ROC for AGE")

F5 = roc(Marked_Recovery ~ Ordered_SEXCD, data=Factors_MR)

plot(F1, main = "ROC for each variables", col="red", type="l", pch = 19)
lines(F2, type="b", col = "royalblue", lty=2, pch = 16)
lines(F3, type = "b", col = "turquoise4", lty =3, pch = 17)
lines(F4, col = "lightcyan4", type ="l", lty=5)
lines(F5, type = "b", col = "black", lty=4, pch = 18)
legend(0.37,0.5, legend=c("Factor 1 (Anemia)", "AIS grades", "Spinal Levels", "Age", "Sex"),
       lty = c(1,2,3,5,4), col = c("red", "royalblue", "turquoise4", "lightcyan4", "black"), 
       pch = c(19,16,17,18), title = "Line Types", cex =0.8, box.lty = 0)
legend(0, 0.5, legend = c("Factor 1 (Anemia) = 63%", "AIS grades = 82%", "Spinal Levels = 60%",
                           "Age = 52%", "Sex = 49%"), title = "Area Under the Curves", cex = 0.8,
       box.lty = 0)

#R^2 for logistics! -> applicable for binomial family! -> change family to binomial
library(pscl) #look for McFadden R^2 scores, ranges from 0-1, closer to 0: no predictive power
pR2(glm_MR)

#Variable importance 
library(caret)
varImp(glm_MR)

#Looking for relative important of each individual variable
#This uses Wald test
library(survey)
regTermTest(glm_MR, "Factor1_w0_flipped")

#Fitted vs residuals plot
qqnorm(resid(glm_MR))
qqline(resid(glm_MR))

#Plots of GLM for odds ratio
library(sjPlot)
sjp.glm(glm_MR)

#Turn estimates in glm output in odds ratio
exp(cbind(Odds=coef(glm_MR), confint(glm_MR)))


plot(glm_MR)

sjp.glm(glm_MR, type = "ma")

#Figures for Marked Recovery 
ggplot(data= subset(Factors_MR, !is.na(Marked_Recovery)), aes(x=Marked_Recovery, y=Factor1_w0))+ 
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Factor 1 (Anemia) Across ASI in Marked Recovery") + 
  facet_grid(SPLVL1 ~ ASIMPC01_A, scales = "free_x")

#DEATH
#Not quite sure if NA means Not Reported or ALIVE!?
# -> Assuming NA means ALIVE
Liver_columns <- c("AST00","ALT00", "CK000")
Liver <- cfa_week0[Liver_columns]

Bio_columns <- c("RBC00", "HGB00", "HCT00", "STP00", "ALB00", "CAB00", "CHO00", "CHC00", "AST00", "ALT00", "CK000",
                 "BUN00","BC900","BUA00","P0400","MCV00","MCH00", "DJN00", "DBN00", "DDN00")
Bio <- subset(cfa_week0, select=(Bio_columns))

#Removing outliers in Factor 3 & 4
Death <- Factors_Overtime_Anemia$DEATHRPD
Death <- as.character(Death)
Death[is.na(Death)] <- 0
Death <- as.factor(Death)
Death <- as.data.frame(Death)
Factors_Overtime_Death <- cbind(Factors_Overtime, Death, Bio, Anemia)
Factors_Overtime_Death <- Factors_Overtime_Death %>% mutate_if(is.character,as.factor)

Factor4_w0_no_outlier <- Factors_Overtime_Death$Factor4_w0
Factor4_w0_no_outlier <- as.data.frame(Factor4_w0_no_outlier)
Factor4_w0_no_outlier <- outlierKD(Factor4_w0_no_outlier, Factor4_w0_no_outlier)

Factor3_w0_no_outlier <- Factors_Overtime_Death$Factor3_w0
Factor3_w0_no_outlier <- as.data.frame(Factor3_w0_no_outlier)
Factor3_w0_no_outlier <- outlierKD(Factor3_w0_no_outlier, Factor3_w0_no_outlier)

#GLM 
glm_death <- glm(Death ~ Factor3_w0 + Factor4_w0 + AGE + SPLVL1 + ASIMPC01_A + SEXCD, data = Factors_Overtime_Death, family = "quasibinomial")
summary(glm_death)

drop1(glm_death, test = "Chi")
stepAIC(glm_death)

glm_death_reduced <- glm(Death ~ Factor3_w0 + Factor4_w0_no_outlier + AGE + ASIMPC01_A, data = Factors_Overtime_Death, family = "quasibinomial")
summary(glm_death_reduced)

anova(glm_death, glm_death_reduced, test="Chisq")

sjp.glm(glm_death_reduced)

Factor3_Death <- ggplot(data= subset(Factors_Overtime_Death, !is.na(ASIMPC01_A)), aes(x=Death, y=Factor3_w0, fill = Death))+ 
  geom_boxplot() + 
  coord_cartesian(ylim=c(-0.85,2)) +
  scale_fill_brewer(name="Death",
                    labels=c("No", "Yes")) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("Factor 3 (Liver) Across ASI grades in Death") + 
  ylab("Factor Liver") + 
  facet_grid(. ~ ASIMPC01_A, scales = "free_x") 

ALT_Death <- ggplot(data= subset(Factors_Overtime_Death, !is.na(ASIMPC01_A)), aes(x=Death, y=ALT00))+ 
  geom_boxplot() +
  coord_cartesian(ylim=c(0,120)) +
  theme_bw() +
  ggtitle("ALT Across ASI grades in Death") + 
  facet_grid(SPLVL1 ~ ASIMPC01_A, scales = "free_x")

AST_Death <- ggplot(data= subset(Factors_Overtime_Death, !is.na(ASIMPC01_A)), aes(x=Death, y=AST00))+ 
  geom_boxplot() + 
  coord_cartesian(ylim=c(0,200)) + 
  theme_bw() +
  ggtitle("AST Across ASI grades in Death") + 
  facet_grid(SPLVL1 ~ ASIMPC01_A, scales = "free_x")

CK_Death <- ggplot(data= subset(Factors_Overtime_Death, !is.na(ASIMPC01_A)), aes(x=Death, y=CK000))+ 
  geom_boxplot() + 
  coord_cartesian(ylim=c(0,6000)) +
  theme_bw() +
  ggtitle("Creatine Phosphokinase Across ASI grades in Death") + 
  facet_grid(SPLVL1 ~ ASIMPC01_A, scales = "free_x")

multiplot(Factor3_Death, ALT_Death, AST_Death, CK_Death, cols = 2)

Factor4_Death <- ggplot(data= subset(Factors_Overtime_Death, !is.na(ASIMPC01_A)), aes(x=Death, y=Factor4_w0, fill=Death))+ 
  geom_boxplot() + 
  coord_cartesian(ylim=c(-1,1.85)) + 
  scale_fill_brewer(name="Death",
                    labels=c("No", "Yes")) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle("Factor 4 (Kidneys) Across ASI grades in Death") + 
  ylab("Factor Kidneys") +
  facet_grid(. ~ ASIMPC01_A, scales = "free_x")

BUN_Death <- ggplot(data= subset(Factors_Overtime_Death, !is.na(ASIMPC01_A)), aes(x=Death, y=BUN00))+ 
  geom_boxplot() + 
  coord_cartesian(ylim=c(0,45)) + 
  theme_bw() +
  ggtitle("Urea Nitrogen Across ASI grades in Death") + 
  facet_grid(SPLVL1 ~ ASIMPC01_A, scales = "free_x")

BC9_Death <- ggplot(data= subset(Factors_Overtime_Death, !is.na(ASIMPC01_A)), aes(x=Death, y=BC900))+ 
  geom_boxplot() + 
  coord_cartesian(ylim=c(0,2)) + 
  theme_bw() +
  ggtitle("Creatinine Across ASI grades in Death") + 
  facet_grid(SPLVL1 ~ ASIMPC01_A, scales = "free_x")

BUA_Death <- ggplot(data= subset(Factors_Overtime_Death, !is.na(ASIMPC01_A)), aes(x=Death, y=BUA00))+ 
  geom_boxplot() + 
  coord_cartesian(ylim=c(0,9)) + 
  theme_bw() +
  ggtitle("Uric Acid Across ASI grades in Death") + 
  facet_grid(SPLVL1 ~ ASIMPC01_A, scales = "free_x")

multiplot(Factor4_Death, BUN_Death, BC9_Death, BUA_Death, cols = 2)

multiplot(Factor3_Death, Factor4_Death)
#RBC drugs 
RBC_drugs$RBC_Drug[is.na(RBC_drugs$RBC_Drug)] <- 0
RBC_drugs$RBC_Drug <- as.factor(RBC_drugs$RBC_Drug)

Other_with_RBC_drug <- merge(Other_Variables, RBC_drugs, by.x=c("PTID"), by.y=c("PTID"))

Factors_Overtime_RBC_drugs <- do.call("cbind", list(MS52_raw, MS00_raw, p_week0, p_week1, p_week2, p_week4, p_week8, p_week52, Other_with_RBC_drug, Anemia, RBC, Death))
Factors_Overtime_RBC_drugs <- Factors_Overtime_RBC_drugs %>% mutate_if(is.character,as.factor)
Factors_Overtime_RBC_drugs$SEXCD <- as.factor(Factors_Overtime_RBC_drugs$SEXCD)

lm_RBC_drugs <- lm(TLAMS52 ~ RBC_Drug, data = Factors_Overtime_RBC_drugs)
lm_factor1 <- lm(TLAMS52 ~ Factor1_w0 + AGE + SPLVL1 + ASIMPC01_A + RBC_Drug + SEXCD, data = Factors_Overtime_RBC_drugs)
lm_reduced_factor1 <- lm(TLAMS52 ~ Factor1_w0 + AGE + SPLVL1 + ASIMPC01_A + RBC_Drug, data = Factors_Overtime_RBC_drugs)
lm_RBC_factor1 <- lm(TLAMS52 ~ Factor1_w0 + RBC_Drug, data = Factors_Overtime_RBC_drugs)
anova(lm_reduced_factor1, lm_factor1, test = "Chisq")

summary(lm_RBC_drugs)
summary(lm_factor1)
sjp.lm(lm_RBC_drugs)

Factors_Overtime_RBC_drugs$RBC_and_Drugs <- ifelse(Factors_Overtime_RBC_drugs$Anemia == 0 & Factors_Overtime_RBC_drugs$RBC_Drug == 0, 0,
                                                   ifelse(Factors_Overtime_RBC_drugs$Anemia == 0 & Factors_Overtime_RBC_drugs$RBC_Drug == 1, 1,
                                                          ifelse(Factors_Overtime_RBC_drugs$Anemia == 1 & Factors_Overtime_RBC_drugs$RBC_Drug == 0, 2, 3)))

Factors_Overtime_RBC_drugs$RBC_and_Drugs <- as.factor(Factors_Overtime_RBC_drugs$RBC_and_Drugs)

step <- stepAIC(lm_factor1, direction = "both")
summary(step)

# Bootstrap Measures of Relative Importance (100 samples) 
library(relaimpo)

calc.relimp(lm_reduced_factor1,type=c("lmg"),
            rela=TRUE)

boot <- boot.relimp(lm_reduced_factor1, b = 100, type = c("lmg"), rank = TRUE, 
                    diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result

plot(booteval.relimp(boot,sort=TRUE)) # plot result

#Relative risk of deaths in Death across RBC drugs
tab <- matrix(c(580,30,170,17),byrow=TRUE,nrow=2)
epitab(tab,method="riskratio")

#Relative risk of deaths in Death across Anemia
tab <- matrix(c(191,13,490,32),byrow=TRUE,nrow=2)
epitab(tab,method="riskratio")

#Lookinng at how Factor1 vs Marked Recovery 
Factor_vs_MR <- ggplot(data= subset(Factors_MR, !is.na(Marked_Recovery)), aes(x=Marked_Recovery, y=Factor1_w0))+ 
  geom_boxplot() + 
  coord_cartesian(ylim=c(-2,3)) +
  theme_bw() +
  ggtitle("Factor 1 (Anemia) Across ASI grades in Marked Recovery") + 
  facet_grid(SPLVL1 ~ ASIMPC01_A, scales = "free_x") 

MS00_vs_MR <- ggplot(data= subset(Factors_MR, !is.na(Marked_Recovery)), aes(x=Marked_Recovery, y=TLAMS00))+ 
  geom_boxplot() +
  theme_bw() +
  facet_grid(. ~ ASIMPC01_A, scales = "free_x") 

MS01_vs_MR <- ggplot(data = subset(Factors_MR, !is.na(Marked_Recovery)),aes(x=Marked_Recovery, y=TLAMS01))+
  geom_boxplot() +
  theme_bw() +
  facet_grid(. ~ ASIMPC01_A, scales = "free_x")

Factor_exclude_MS5 <- ggplot(data = subset(Factors_MR, TLAMS00 <= 5 &!is.na(Marked_Recovery)),aes(x=Marked_Recovery, y=Factor1_w0))+
  geom_boxplot() +
  coord_cartesian(ylim=c(-2,3)) + 
  theme_bw() +
  ggtitle("Factor 1 (Anemia) Across ASI grades in Marked Recovery", subtitle = "(Excluding Motor Score > 5)") +
  facet_grid(SPLVL1 ~ ASIMPC01_A, scales = "free_x")
  
multiplot(Factor_vs_MR, Factor_exclude_MS5)

#Comparing C's between Marked Recovery of Factor 1
MR <- c("Factor1_w0", "Marked_Recovery", "ASIMPC01_A", "TLAMS00")
Factors_MR_Cs <- Factors_MR[MR]
Factors_MR_Cs <- subset(Factors_MR_Cs, ASIMPC01_A == "C")
t.test(Factor1_w0~Marked_Recovery, data = Factors_MR_Cs)

ggplot(data = subset(Factors_Overtime_RBC_drugs, !is.na(RBC_and_Drugs)), aes(x = RBC_and_Drugs, y = TLAMS52, fill = RBC_and_Drugs)) +
 geom_boxplot() +
 theme_bw() +
 theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
 ggtitle("RBC Drugs and Factor 1 (Anemia Factor)") +
 scale_fill_manual(name="Anemia \n and RBC Drugs",
                      breaks=c("0", "1", "2", "3"),
                      labels=c("Non-Anemic and Not taking Drugs", 
                               "Non-Anemic and Taking Drugs", 
                               "Anemic and Not taking Drugs",
                               "Anemic and Taking Drugs"),
                     values = c("lightgoldenrod1", "darkgoldenrod", "paleturquoise", "paleturquoise4")) + 
  facet_grid(. ~ Marked_Recovery, scales = "free_x")

ggplot(data = Factors_Overtime_RBC_drugs) + 
  geom_mosaic(aes(x = product(RBC_and_Drugs, Marked_Recovery), fill=factor(RBC_and_Drugs)), na.rm=TRUE, offset=0.015) +
  theme_bw() +
  ggtitle("RBC Drugs and MR") +
  scale_fill_manual(name="Anemia \nand RBC Drugs",
                    breaks=c("0", "1", "2", "3"),
                    labels=c("Non-Anemic and Not taking Drugs", 
                             "Non-Anemic and Taking Drugs", 
                             "Anemic and Not taking Drugs",
                             "Anemic and Taking Drugs"),
                    values = c("lightgoldenrod1", "darkgoldenrod", "paleturquoise", "paleturquoise4")) + 
  facet_grid(. ~ Marked_Recovery, scales = "free_x")

#if you are anemic and you get RBC, do you do better than if you are anemic and do not get RBC?
Anemia_RBC_only <- subset(Factors_Overtime_RBC_drugs, RBC_and_Drugs %in% c("2","3"))

glm_drug <- glm(Marked_Recovery ~ RBC_and_Drugs, data = Anemia_RBC_only, family = quasibinomial)
summary(glm_drug)

sjp.glm(glm_drug)

#Survival Analysis for Death 
#Interpretation: positive coef means that risk of death is higher, and negative means risk is lower
#Example: Coef for Age is 0.043, means older age has higher risk of death
#Coef for Factor 3_cutoff is -0.57 means that higher liver enzymes = lower risk of death 
#Hazard Ratio = exponentiated coef
#HR for Factor 3_cutoff is 0.56 means that higher liver enzymes reduces hazard ratio by 0.56
#In other words, HR for Factor3_cutoff reduces hazard ratio by 44%


#The survival probability is the probability that an individual survives from the time origin to a specified future time t.
#The hazard is the probability that an individual who is under observation at a time t has an event at that time.

#Survivor function focuses on NOT having an event & the hazard function focuses on the event occurring.
Other_Death$Survival_days[is.na(Other_Death$Survival_days)] <- 365

Factors_Overtime_Death_Sur <- merge(Factors_Overtime_Death, Other_Death, by.x=c("PTID"), by.y=c("ptid"))

Factors_Overtime_Death_Sur$Survival_Rates <- ifelse(Factors_Overtime_Death_Sur$Death == 1, 0,1)
Factors_Overtime_Death_Sur$Survival_Rates <- as.factor(Factors_Overtime_Death_Sur$Survival_Rates)

Factors_Overtime_Death_Sur$Factor3_w0_cutoff <- ifelse(Factors_Overtime_Death_Sur$Factor3_w0_flipped <= 0.74253,0,1)
Factors_Overtime_Death_Sur$Factor3_w0_cutoff <- as.factor(Factors_Overtime_Death_Sur$Factor3_w0_cutoff)

Factors_Overtime_Death_Sur$Factor4_w0_cutoff <- ifelse(Factors_Overtime_Death_Sur$Factor4_w0_flipped <= 0.91,0,1)
Factors_Overtime_Death_Sur$Factor4_w0_cutoff <- as.factor(Factors_Overtime_Death_Sur$Factor4_w0_cutoff)

Factors_Overtime_Death_Sur$Severity <- ifelse(Factors_Overtime_Death_Sur$ASIMPC01_A == "A", 1,2)
Factors_Overtime_Death_Sur$Severity <- as.factor(Factors_Overtime_Death_Sur$Severity)  

Factors_Overtime_Death_Sur$SPLVL1_numeric <- ifelse(Factors_Overtime_Death_Sur$SPLVL1 == "C", 1,2)
Factors_Overtime_Death_Sur$SPLVL1_numeric <- as.factor(Factors_Overtime_Death_Sur$SPLVL1_numeric)  

Factors_Overtime_Death_Sur$Factor3_w0_flipped <- Factors_Overtime_Death_Sur$Factor3_w0+1
Factors_Overtime_Death_Sur$Factor4_w0_flipped <- Factors_Overtime_Death_Sur$Factor4_w0+1

library(survival)

Survival_Factor3 <- survfit(Surv(Survival_days, Survival_Rates==0) ~ Factor3_w0_cutoff, data = Factors_Overtime_Death_Sur)
summary(Survival_Factor3)$table

Survival_Factor4 <- survfit(Surv(Survival_days, Survival_Rates == 0) ~ Factor4_w0_cutoff, data = Factors_Overtime_Death_Sur)
summary(Survival_Factor4)$table

library(survminer)
ggsurvplot(Survival_Factor3,      
  data = Factors_Overtime_Death_Sur,  
  risk.table = "abs_pct",
  pval = TRUE,             
  conf.int = TRUE, 
  xlim = c(0,500),     
  risk.table.y.text.col = T, 
  risk.table.y.text = FALSE,
  ncensor.plot = TRUE,
  surv.median.line = "hv",
  loglog   = FALSE,                 
  logt     = FALSE) 

ggsurvplot(Survival_Factor4,      
           data = Factors_Overtime_Death_Sur,  
           risk.table = "abs_pct",
           pval = TRUE,             
           conf.int = TRUE, 
           xlim = c(0,500),     
           risk.table.y.text.col = T, 
           risk.table.y.text = FALSE,
           ncensor.plot = TRUE,
           surv.median.line = "hv",
           loglog   = FALSE,                 
           logt     = FALSE) 

survdiff(Surv(Survival_days,Survival_Rates == 0)~Factor3_w0_cutoff, data = Factors_Overtime_Death_Sur)

#Univariate models
covariates <- c("Factor1_w0","Factor2_w0","Factor3_w0", "Factor4_w0","Factor5_w0","Factor6_w0")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(Survival_days, Survival_Rates == 0)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = Factors_Overtime_Death_Sur)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         return(exp(cbind(coef(x),confint(x))))
                       })

Survival_Factor3_all <- as.data.frame(univ_results, check.names = TRUE)

write.csv(Survival_Factor3_all, 'Univariate_Cox.csv')
write.table(Survival_Factor3_all, 'Univariate_Cox.txt')

#Multivariate 
Hazard_Factor3 <- coxph(Surv(Survival_days, Survival_Rates==0)~ Factor3_w0_flipped + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Death_Sur)
summary(Hazard_Factor3)

Hazard_Factor4 <- coxph(Surv(Survival_days, Survival_Rates==0)~ Factor4_w0_flipped + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Death_Sur)
summary(Hazard_Factor4)

Hazard_Factors <- coxph(Surv(Survival_days, Survival_Rates==0)~ Factor3_w0_flipped + Factor4_w0_flipped + AGE + SEXCD + SPLVL1 + ASIMPC01_A, dist = "weibull", data = Factors_Overtime_Death_Sur)
summary(Hazard_Factors)

library(survminer)
ggforest(Hazard_Factor3, data = Factors_Overtime_Death_Sur)
library(forestmodel)
forest_model(Hazard_Factor3)

ggsurvplot(survfit(Hazard_Factor3), palette = "#2E9FDF",
           ggtheme = theme_minimal(), data = Factors_Overtime_Death_Sur)

#HR curves for continous predictors
library(smoothHR)
hr_plot <- smoothHR(data = Factors_Overtime_Death_Sur, coxfit = Hazard_Factors)
plot(hr_plot, predictor = "Factor3_w0_flipped", prob = 0, conf.level = 0.95)

#Diagnostics for Cox Proportinal Hazard

#Test the proportional hazards assumption for each covariate included by scaled Schoenfeld residuals
#In principle, the Schoenfeld residuals are independent of time. 
#The proportional hazard assumption is supported by a non-significant relationship between residuals and time, 
#and refuted by a significant relationship.
Cox_test <- cox.zph(Hazard_Factor3) 
#A plot that shows a non-random pattern against time is evidence of violation of the PH assumption
ggcoxzph(Cox_test)

#testing influential observations
#The deviance residual is a normalized transform of the martingale residual. 
#These residuals should be roughly symmetrically distributed about zero with a standard deviation of 1

ggcoxdiagnostics(Hazard_Factor3, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

#Positive values correspond to individuals that "died too soon" compared to expected survival times.
#Negative values correspond to individual that "lived too long".
#Very large or small values are outliers, which are poorly predicted by the model.
ggcoxdiagnostics(Hazard_Factor3, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())

#testing for linearity 
#Plotting the Martingale residuals against continuous covariates to detect nonlinearity
#patterns in the plot may suggest that the variable is not properly fit
# Fitted lines with lowess function should be linear to satisfy the Cox proportional hazards model assumptions.
#Factor 4 is not linear!!!!

Martingale_resids <- residuals(Hazard_Factor3, type = c("martingale"), collapse = Factors_Overtime_Death_Sur$Factor3_w0)
Martingale_resids <- as.data.frame(Martingale_resids)
Martingale_resids$new_Id <- 1:784

Factor3_w0 <- Factors_Overtime_Death_Sur$Factor3_w0
Factor3_w0 <- as.data.frame(Factor3_w0)
Factor3_w0$new_Id <- 1:797

Factor4_w0 <- Factors_Overtime_Death_Sur$Factor4_w0
Factor4_w0 <- as.data.frame(Factor4_w0)
Factor4_w0$new_Id <- 1:797

Age <- Factors_Overtime_Death_Sur$AGE
Age <- as.data.frame(Age)
Age$new_Id <- 1:797

Resids_Factors <- Reduce(function(x, y) merge(x, y, all=TRUE, by="new_Id"), 
                 list(Factor3_w0, Factor4_w0, Martingale_resids, Age))

ggplot(data = Resids_Factors, aes(x=Factor4_w0, y=Martingale_resids))+
  geom_point() +
  geom_smooth(method= "loess")


summary(gam(Survival_Rates ~ Factor4_w0 + AGE + ASIMPC01_A + SEXCD + SPLVL1, family = cox.ph(), data = Factors_Overtime_Death_Sur))

#Check for each individual variables for the Cox Ph
aa_fit <-aareg(Surv(Survival_days, Survival_Rates == 0) ~ Factor3_w0 + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Death_Sur)
autoplot(aa_fit)

#Weibull Distribution - Accelerated Failure Time model (AFT)
library(SurvRegCensCov)
Wei_Factor3 <- survreg(Surv(Survival_days, Survival_Rates == 0) ~ Factor3_w0_flipped + Factor4_w0_flipped + AGE + SEXCD + SPLVL1 + ASIMPC01_A, dist = "weibull", data = Factors_Overtime_Death_Sur)
Converted_Wei <- ConvertWeibull(Wei_Factor3, conf.level = 0.95)

Wei_Factor<- WeibullReg(Surv(Survival_days, Survival_Rates == 0) ~ Factor3_w0 + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Death_Sur) #Same as above

#Weilbull Diagnostics 
WeibullDiag(Surv(Survival_days, Survival_Rates == 0) ~ SEXCD, data = Factors_Overtime_Death_Sur)

library(eha)
Weibull_Factor3_alt <- weibreg(Surv(Survival_days, Survival_Rates == 0) ~ ALT00 + AGE + SEXCD + ASIMPC01_A + SPLVL1, data = Factors_Overtime_Death_Sur)

plot(Weibull_Factor3_alt,fn=c("sur"),new.data = c(17,35,2,1,1))
plot(Weibull_Factor3_alt,fn=c("sur"),new.data = c(70,35,2,1,1))

#Goodness-of-fit test
phreg.factor3 <- phreg(Surv(Survival_days, Survival_Rates == 0) ~ Factor3_w0_flipped + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Death_Sur, dist = "weibull")
coxreg.factor3 <- coxreg(Surv(Survival_days, Survival_Rates == 0) ~ Factor3_w0_flipped + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Death_Sur)
check.dist(coxreg.factor3, phreg.factor3)

#Importance of covariates 
psm.factor3 <- psm(Surv(Survival_days, Survival_Rates == 0) ~ Factor3_w0 + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Death_Sur)

#TABLES
sjt.lm(lm_week0_factor1, lm_week0_factor2, lm_week0_factor3,
       p.numeric = FALSE,
       depvar.labels = c("Factor1 (Anemia)", "Factor2 (Liver)","Factor3 (Kidneys)"),
       string.est = "Estimate", 
       string.ci = "C.I",  
       show.header = TRUE, 
       sep.column = FALSE,
       string.dv = "Multivariate Analysis (Adjusted) for Motor Scores with Different Factors", 
       string.pred = "Coefficients",
       cell.spacing = 0.1,
       pred.labels = c("Factor1 (Anemia)", "Sex (Male)", "Spinal level (Throacic)", "Age",
                       "AIS-A vs B", "AIS-A vs C", "AIS-A vs D", "Factor2 (Nutrition)", "Factor3 (Liver)"),
       group.pred = FALSE)

sjt.lm(lm_week0_factor1, lm_week0_factor3, lm_week0_factor4, lm_week0_factor5, lm_week0_factor6, lm_week0_factor2,
        p.numeric = FALSE,
        depvar.labels = c("Factor Anemia", "Factor Liver",  "Factor Kidneys", 
                          "Factor RBC Indices", "Factor WBC", "Factor Others"),
        string.est = "Estimate",
        string.ci = "C.I",  
        show.header = TRUE, 
        sep.column = FALSE,
        string.dv = "Univariate Analysis for LEMS with Different Factors", 
        string.pred = "Coefficients",
        cell.spacing = 0.1,
        pred.labels = c("Factor Anemia", "Factor Liver", "Factor Kidneys", 
                        "Factor RBC Indices", "Factor WBC", "Factor Others"),
        group.pred = FALSE)

sjt.lm(lm_week0_factor1, lm_week0_factor1_updated,
       p.numeric = TRUE,
       depvar.labels = c("Initial Model", "Updated Model"),
       string.est = "Estimate", 
       string.ci = "C.I",  
       show.header = TRUE, 
       sep.column = FALSE,
       string.dv = "Multivariate Analysis (Adjusted) Between 2 Models", 
       string.pred = "Coefficients",
       cell.spacing = 0.1,
       group.pred = FALSE,
       pred.labels = c("Factor1 (Anemia)", "Sex (Male)", "Spinal level (Throacic)", "Age",
                       "AIS-A vs B", "AIS-A vs C", "AIS-A vs D"))

sjt.glm(glm_MR_factor1,  glm_MR_factor3, glm_MR_factor4, glm_MR_factor5, glm_MR_factor6, glm_MR_factor2,
        p.numeric = FALSE,
        depvar.labels = c("Factor Anemia", "Factor Liver", "Factor Kidneys",
                           "Factor RBC Indices", "Factor WBC", "Factor Others"),
        exp.coef = TRUE,
        string.ci = "C.I",  
        show.header = TRUE, 
        sep.column = FALSE,
        string.dv = "Univariate Analysis for Marked Recovery with Different Factors", 
        string.pred = "Coefficients",
        cell.spacing = 0.1,
        pred.labels = c("Factor Anemia", "Factor Liver", "Factor Kidneys",
                        "Factor RBC Indices", "Factor WBC", "Factor Others"),
        group.pred = FALSE)

#Palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(data = subset(Factors_MR, !is.na(Marked_Recovery)), aes(x=Marked_Recovery, y=Factor1_w0_flipped, fill = Marked_Recovery))+
  geom_boxplot() +
  scale_fill_brewer(name="Marked Recovery",
                    labels=c("No", "Yes")) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("Factor 1 (Anemia) and Marked Recovery") +
  ylab("Factor 1 (Anemia)") +
  facet_grid(. ~ ASIMPC01_A, scales = "free_x")

lm_week0_factor1_A <- lm(TLAMS52 ~  Factor1_w0 + SPLVL1 + AGE, ASIMPC01_A == "A", data = Factors_W0)
lm_week0_factor1_B <- lm(TLAMS52 ~  Factor1_w0 + SPLVL1 + AGE, ASIMPC01_A == "B", data = Factors_W0)
lm_week0_factor1_c <- lm(TLAMS52 ~  Factor1_w0 + SPLVL1 + AGE, ASIMPC01_A == "C", data = Factors_W0)
lm_week0_factor1_D <- lm(TLAMS52 ~  Factor1_w0 + SPLVL1 + AGE, ASIMPC01_A == "D", data = Factors_W0)

sjt.lm(lm_week0_factor1_A, lm_week0_factor1_B, lm_week0_factor1_c, lm_week0_factor1_D,
       p.numeric = FALSE,
       depvar.labels = c("A only", "B only","C only","D only"),
       string.est = "Estimate", 
       string.ci = "C.I",  
       show.header = TRUE, 
       sep.column = FALSE,
       string.dv = "Multivariate Analysis (Adjusted) Stratified by AIS grades", 
       string.pred = "Coefficients",
       cell.spacing = 0.1,
       group.pred = FALSE,
       pred.labels = c("Factor1 (Anemia)", "Age", "Spinal level (Throacic)", "Age",
                       "AIS-A vs B", "AIS-A vs C", "AIS-A vs D"))

ggplot(data = subset(Factors_MR, !is.na(Marked_Recovery)), aes(x=Marked_Recovery, y=Factor2_w0))+
  geom_boxplot() +
  ggtitle("Factor 2 (Nutrition) and Marked Recovery") +
  ylab("Factor1 (Anemia)") +
  theme_bw()

sjt.glm(glm_MR, glm_MR_2,
        p.numeric = FALSE,
        depvar.labels = c("Anemia Model (Adjusted)", "Nutrition Model (Adjusted)"),
        exp.coef = TRUE,
        string.ci = "C.I",  
        show.header = TRUE, 
        sep.column = FALSE,
        string.dv = "Multivariate Analysis (Adjusted) for Marked Recovery with Different Factors", 
        string.pred = "Coefficients",
        cell.spacing = 0.1,
        pred.labels = c("Factor1 (Anemia)", "Sex (Male)", "Spinal level (Throacic)", "Age",
                        "AIS-A vs B", "AIS-A vs C", "AIS-A vs D", "Factor2 (Nutrition)"),
        group.pred = FALSE)

sjt.glm(glm_MR, glm_MR_reduced,
        p.numeric = TRUE,
        depvar.labels = c("Initial Model (Adjusted)", "Updated Model (Adjusted)"),
        exp.coef = TRUE,
        string.ci = "C.I",  
        show.header = TRUE, 
        sep.column = FALSE,
        string.dv = "Multivariate Analysis (Adjusted) for Marked Recovery with Different Factors", 
        string.pred = "Coefficients",
        cell.spacing = 0.1,
        pred.labels = c("Factor1 (Anemia)", "Sex (Male)", "Spinal level (Throacic)", "Age",
                        "AIS-A vs B", "AIS-A vs C", "AIS-A vs D", "Factor2 (Nutrition)"),
        group.pred = FALSE)


#WHY LOWER LIVER ENZYMES = DEATH?
cfa2_week0_all <- cbind(Factors_Overtime_Death, cfa_week0)
cfa2_week0_all$Factor3_w0_flipped <- cfa2_week0_all$Factor3_w0+1
cfa2_week0_all$Factor4_w0_flipped <- cfa2_week0_all$Factor4_w0+1
cfa2_week0_all$Liver_enzymes <- ifelse(cfa2_week0_all$Factor3_w0_flipped <= 0.74, 1, 0)
cfa2_week0_all$Liver_enzymes <- as.factor(cfa2_week0_all$Liver_enzymes)

ggplot(data = cfa2_week0_all, aes(x=MCH00, y=ALT00))+
  geom_point()

cfa2_week0_all_dead <- subset(cfa2_week0_all, Death=="1")

cfa2_week0_all_alive <- subset(cfa2_week0_all, Death=="0")

cfa2_week0_all_liver <- subset(cfa2_week0_all, Liver_enzymes=="1")

#Tables and graphs 
forest_model(lm_week0_Anemia)
forest_model(glm_MR)
forest_model(Hazard_Factors)

#ETR tables
ETR <- Converted_Wei$ETR
ETR <- as.data.frame(ETR)
ETR$variable <- c("Factor3", "Factor4", "Age", "Sex", "Spinal Level Injury Thoracic", "AIS - B", "AIS - C", "AIS - D")


ggplot(ETR, aes(x=variable, y=ETR, ymin=LB, ymax=UB))+
  geom_pointrange()+
  geom_hline(yintercept = 1, linetype=2)+
  coord_flip()+
  xlab('Variable')

ETR$se <- log(ETR$UB)-log(ETR$LB)/(2*1.96)
res  <- rma(yi=ETR$ETR, sei=ETR$se, method="FE")
forest(res)

forest_rma(res)
row_names <- cbind(c("Subgroup","\n",ETR$variable), 
                                c("4-Yr Cum. Event Rate\n PCI","\n",ETR$ETR))

forestplot(labeltext = c("Factor Liver", "Factor Kidney", "Age", "Sex", "Spinal Level Injury", "AIS-B", 
                         "AIS-c", "AIS-D"), 
          mean=ETR$ETR, 
          lower=ETR$LB, upper=ETR$UB,
          title="Hazard Ratio", 
          lwd.ci = 2)

#Heat plot
cfa_week0_variables <- cfa_week0[c(-1,-2,-3,-40)]


Cor_cfa <- cor(cfa_week0_variables, use="pairwise.complete.obs")
Cor_cfa <- as.matrix(Cor_cfa)

corrplot(Cor_cfa, method = "color", type = "upper", order = "FPC",
         addCoef.col = "black", number.cex=0.35)
