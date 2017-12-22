
#Scaling data of week 1
cfa_week1 <- cbind(MS52_raw, cfa_week1)
cfa1_week1 <- scale(cfa_week1, center = TRUE, scale = TRUE)
cfa1_week1 <- data.frame(cfa1_week1)
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
cfa_week2 <- cbind(MS52_raw, cfa_week2)
cfa1_week2 <- scale(cfa_week2, center = TRUE, scale = TRUE)
cfa1_week2 <- data.frame(cfa1_week2)
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


fa <- fa(cfa1_week4, nfactors=6, rotate="varimax", fm="ml")
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
cfa1_week8

cfa_week8 <- cbind(MS52_raw, cfa_week8)

cor(cfa_week8, use = "pairwise.complete.obs")

fa <- fa(cfa1_week8[-1], nfactors=6, rotate="varimax", fm="ml")
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
cfa_week52 <- cbind(MS52_raw, cfa_week52)

cfa1_week52 <- scale(cfa_week52, center = TRUE, scale = TRUE)
cfa1_week52 <- data.frame(cfa1_week52)
cfa1_week52

cor(cfa_week52, use="pairwise.complete.obs")

fa <- fa(cfa1_week52[-1], nfactors=6, rotate="varimax", fm="ml")
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


#Raw and Scaled scores for MS52 and MS00
MS52_scaled <- as.data.frame(scale(cfa_week0[1], scale = TRUE, center = TRUE))
MS52_raw <- cfa_week0[1]

MS00_raw <- cfa_week0[2]

#Combined dataframe
Factors_Overtime <- do.call("cbind", list(MS00_raw,MS52_scaled, p_week0, p_week1, p_week2, p_week4, p_week8, p_week52, Other_Variables))

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

UCP_week1 <- ctree(TLAMS52 ~ Factor1_w1 + Factor2_w1 +Factor3_w1 + Factor4_w1 + Factor5_w1 + Factor6_w1 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
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

lm_week0_factor1 <- lm(TLAMS52 ~  Factor1_w0 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
lm_week0_factor2 <- lm(TLAMS52 ~  Factor2_w0 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
lm_week0_factor3 <- lm(TLAMS52 ~  Factor3_w0 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)

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

#Model diagnostics
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

lm_week0_Anemia <- lm(TLAMS52 ~ Factor1_w0 + AGE + SPLVL1 + ASIMPC01_A + SEXCD, data = Factors_Overtime_Anemia)
summary(lm_week0_Anemia)
sjp.lm(lm_week0_Anemia)
drop1(lm_week0_Anemia, test = "Chi")

lm_week0_Anemia_inj <- lm(TLAMS52 ~ SEXCD + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Anemia)
summary(lm_week0_Anemia_inj)

lm_week0_Anemia_alone <- lm(TLAMS52 ~ Anemia, data = Factors_Overtime_Anemia)
summary(lm_week0_Anemia_alone)

lm_week0_Anemia_Factor <- lm(TLAMS52 ~ Factor1_w0 + AGE + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Anemia)
summary(lm_week0_Anemia_Factor)


AIC(lm_week0_Anemia)
AIC(lm_week0_Anemia_inj)
AIC(lm_week0_Anemia_alone)
AIC(lm_week0_Anemia_Factor)

anova(lm_week0_Anemia, lm_week0_Anemia_Factor, test = "Chisq")

#Diagnostic
qqnorm(resid(lm_week0_Anemia_Factor))
qqline(resid(lm_week0_Anemia_Factor))

#Diagnostic plot -> looks like we have heteroscedascity and non-normality problem
autoplot(lm_week0_Anemia_Factor, which = 1:6, colour = 'ASIMPC01_A')
plot(lm_week0_Anemia_Factor, add.smooth = TRUE)

#Skewed data
library(moments)
skewness(lm_week0_Anemia_Factor$residuals)

plot(Effect(c("Factor1_w0", "ASIMPC01_A"), lm_week0_Anemia_Factor), multiline = TRUE)

#Logarithm is clumsy and doesn't work
Factors_Overtime_Anemia$log.MS52 <- log(Factors_Overtime_Anemia$TLAMS52 + 1)
Factors_Overtime_Anemia$log.Factor1_w0 <- log(Factors_Overtime_Anemia$Factor1_w0 + 1)

#Graph of lm and diagnostic 
library(sjPlot)
sjp.lm(lm_week0_Anemia_Factor)

#Model Diagnostic as follows: Outliers -> Linearity -> Heteroscedasticiy -> Non-normality

#Outliers
outlierTest(lm_week0_Anemia_Factor)
cooksd <- cooks.distance(lm_week0_Anemia_Factor)
plot(cooksd, pch="*", cex=1, main="Influential Obs by Cooks distance")
abline(h = 4*mean(cooksd, na.rm=T), col="red") 
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")

#Best Outlier test!!! I think - no outliers
source("http://goo.gl/UUyEzD")
b <- outlierKD(Factors_Overtime_Anemia, Factor1_w0)

b <- cbind(Factors_Overtime_Anemia, b)
#No outliers worth removing 

#Check for linearity - Yes
crPlot(lm_week0_Anemia_Factor, "Factor1_w0")

spreadLevelPlot(lm_week0_Anemia_Factor)

#Multicolinearity - No
library(car)
vif(lm_week0_Anemia_Factor)

#Other tests
library(gvlma)
gvmodel <- gvlma(lm_week0_Anemia_Factor) 
summary(gvmodel)

#Normal distribution? -> not
shapiro.test(resid(lm_week0_Anemia_Factor))

#Heteroscasdicity -> yes!
ncvTest(lm_week0_Anemia_Factor)
bptest(lm_week0_Anemia_Factor)


# We can also try BoxCox but with our dataset, boxcox may not be appropriate because of extreme values (0)

# Option A for BoxCox transformation
boxcoxMS <- caret::BoxCoxTrans(Factors_Overtime_Anemia$TLAMS52+1)
print(boxcoxMS)
new_boxcox <- cbind(Factors_Overtime_Anemia, MS52_boxcox=predict(boxcoxMS, Factors_Overtime_Anemia$TLAMS52+1))
lm_new_boxcox <- lm(MS52_boxcox ~  Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = new_boxcox)
summary(lm_new_boxcox)
ncvTest(lm_new_boxcox)
bptest(lm_new_boxcox)
shapiro.test(resid(lm_new_boxcox))
qqnorm(resid(lm_new_boxcox))
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

# => Both options did NOT transform data!!!!!!!!
# ncvTest and bptest are different because bptest is studendized and ncvTest is not
# => which means bptest is more robust than the original one (ncvTest)
# => go with bptest -> null is: homoscedasicity 

#Use this to treat for heteroscadasicity 
library(lmtest)
coeftest(lm_week0_Anemia_Factor,vcov=hccm(lm_week0_Anemia_Factor))
HC <- hccm(lm_week0_Anemia_Factor)

homo <- coeftest(lm_week0_Anemia_Factor, vcov. = HC)
waldtest(lm_week0_Anemia_Factor, vcov = HC)

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

Factors_Overtime_Anemia <- cbind(Factors_Overtime_Anemia, RBC)

UCP_RBC <- partykit::ctree(TLAMS52 ~ RBC00 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
plot(UCP_RBC)

UCP_HGB <- partykit::ctree(TLAMS52 ~ HGB00 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
plot(UCP_HGB)

UCP_HCT <- partykit::ctree(TLAMS52 ~ HCT00 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
plot(UCP_HCT)

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


Factors_MR <- do.call("cbind", list(MS52_raw, p_week0, p_week1, p_week2, p_week4, p_week8, p_week52, Other_Variables, RBC))

# GLM for marked Recovery 
Factors_MR <- Factors_MR %>% mutate_if(is.character,as.factor)
Factors_MR$SEXCD <- as.factor(Factors_MR$SEXCD) 
Factors_MR$AGE <- as.numeric(Factors_MR$AGE)
Factors_MR$Marked_Recovery <- as.factor(Factors_MR$Marked_Recovery)


glm_MR_factor1 <- glm(Marked_Recovery ~ Factor1_w0, data = Factors_MR, family = quasibinomial)
summary(glm_MR_factor1)

glm_MR_inj <- glm(Marked_Recovery ~ SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_MR, family = quasibinomial)
summary(glm_MR_inj)

anova(glm_MR_factor1, glm_MR_reduced, test = "Chisq")

glm_MR <- glm(Marked_Recovery ~ Factor1_w0 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_MR, family = quasibinomial)
summary(glm_MR)

drop1(glm_MR_inj, test = "Chi")

glm_MR_reduced <- glm(Marked_Recovery ~ Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_MR, family = quasibinomial)
summary(glm_MR_reduced)


anova(glm_MR_reduced, glm_MR, test = "Chisq")

#test for dispersion 
P__disp(glm_MR_reduced)

#Ctree
Factors_MR <- subset(Factors_MR, !is.na(Marked_Recovery))

UCP_MR <- ctree(Marked_Recovery ~ Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_MR)
plot(UCP_MR)

#Model diagnostics
crPlots(glm_MR)
vif(glm_MR)

library(boot)
glm.diag.plots (glm_MR_reduced)


library(pROC)
AUC = roc(Marked_Recovery ~ Factor1_w0, data=Factors_MR, plot = TRUE, percent = TRUE, ci = TRUE) 


library(pscl)
pR2(glm_death_reduced)

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


#DEATH
#Not quite sure if NA means Not Reported or ALIVE!?
# -> Assuming NA means ALIVE

Death <- Factors_Overtime_Anemia$DEATHRPD
Death[is.na(Death)] <- 0
Death[Death == 2] <- 1
Death <- as.factor(Death$Death)
Death <- as.data.frame(Death)
Factors_Overtime_Death <- cbind(Factors_Overtime, Death)
Factors_Overtime_Death <- Factors_Overtime_Death %>% mutate_if(is.character,as.factor)

glm_death <- glm(Death ~ Factor1_w0 + AGE + SPLVL1 + ASIMPC01_A + SEXCD, data = Factors_Overtime_Death, family = "quasibinomial")
summary(glm_death)

drop1(glm_death, test = "Chi")
stepAIC(glm_death)

glm_death_reduced <- glm(Death ~ AGE + ASIMPC01_A, data = Factors_Overtime_Death, family = "quasibinomial")
summary(glm_death_reduced)

anova(glm_death, glm_death_reduced, test="Chisq")

sjp.glm(glm_death)

