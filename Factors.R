#Raw and Scaled scores for MS52 and MS00
MS52_scaled <- as.data.frame(scale(cfa_week0[1], scale = TRUE, center = TRUE))
MS52_raw <- cfa_week0[1]

MS00_raw <- cfa_week0[2]

MS01_raw <- cfa_week0[3]

MS04_raw <- MS_only["TLAMS04"]

MS08_raw <- MS_only["TLAMS08"]

MS16_raw <- MS_only["TLAMS16"]

MS26_raw <- MS_only["TLAMS26"]

cfa_week0_only <- cfa_week0[4:39]
cfa_week0_only_scaled <- as.data.frame(scale(cfa_week0[4:39], scale = TRUE, center = TRUE))

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

write.xlsx(Factors_Overtime, 'Factors_Overtime.xlsx')
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
Factors_W0_noMS <- subset(Factors_W0, !is.na(TLAMS52))

#UCP for LMS52
library(party)
UCP_Factors_W0 <- ctree(TLAMS52 ~ Factor1_w0 + AGE + SPLVL1 + ASIMPC01_A, data = Factors_W0_noMS)
plot(UCP_Factors_W0)


#UCP for week 4
Factors_W4 <- do.call("cbind", list(MS52_raw, p_week4, Other_Variables, Factor1))

Factors_W4 <- Factors_W4 %>% mutate_if(is.character,as.factor) #They only take factor values
Factors_W4$SEXCD <- as.factor(Factors_W4$SEXCD) 

Factors_W4_noMS <- subset(Factors_W4, !is.na(TLAMS52))

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
UCP_week0 <- ctree(Marked_Recovery ~ SPLVL1 + SEXCD + AGE + ASIMPC01_A + slope_Factor, data = subset(Factors_Overtime_coefs_update, !is.na(Marked_Recovery)))
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

Factors_Overtime_Anemia <- cbind(Factors_Overtime, Anemia, cfa_week0_only)
Factors_Overtime_Anemia_scaled <- cbind(Factors_Overtime, Anemia, cfa_week0_only_scaled)

Factors_Overtime_Anemia <- subset(Factors_Overtime_Anemia, !is.na(ASIMPC01_A))

Factors_Overtime_Anemia <- subset(Factors_Overtime_Anemia, !is.na(Anemia))
Factors_Overtime_Anemia <- subset(Factors_Overtime_Anemia, !is.na(TLAMS52))

Factors_Overtime_Anemia <- Factors_Overtime_Anemia %>% mutate_if(is.character,as.factor)

Factors_Overtime_Anemia$SEXCD <- as.factor(Factors_Overtime_Anemia$SEXCD) 
Factors_Overtime_Anemia$AGE <- as.numeric(Factors_Overtime_Anemia$AGE)
Factors_Overtime_Anemia$DEATHRPD <- as.factor(Factors_Overtime_Anemia$DEATHRPD)
Factors_Overtime_Anemia$Factor1_w0_flipped <- Factors_Overtime_Anemia$Factor1_w0+1
Factors_Overtime_Anemia$Severity <- ifelse(Factors_Overtime_Anemia$ASIMPC01_A == "A", 1,
                                      ifelse(Factors_Overtime_Anemia$ASIMPC01_A == "B", 1,
                                             ifelse(Factors_Overtime_Anemia$ASIMPC01_A == "C", 0, 0)))
Factors_Overtime_Anemia$Severity <- as.factor(Factors_Overtime_Anemia$Severity)

Factors_Overtime_Anemia_surg <- merge(surg_vars, Factors_Overtime_Anemia, by.x=c("ptid"), by.y=c("PTID"))
Factors_Overtime_AIS <- merge(Factors_Overtime_Anemia_surg, new_AIS, by.x=c("ptid"), by.y=c("PTID"))
Factors_Overtime_AIS$AIS_4 <- as.factor(Factors_Overtime_AIS$AIS_4)
Factors_Overtime_AIS$Severity_w4 <- ifelse(Factors_Overtime_AIS$AIS_4 == "A", 1,
                                           ifelse(Factors_Overtime_AIS$AIS_4 == "B", 1,
                                                  ifelse(Factors_Overtime_AIS$AIS_4 == "C", 0, 0)))



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

lm_week0_Anemia_variables <- lm(TLAMS52 ~ RBC00*HGB00*HCT00 + AGE + SEXCD  + ASIMPC01_A + SPLVL1, data = Factors_Overtime_Anemia_scaled)
summary(lm_week0_Anemia_variables)

#Predictive Model 
#Split data in half to build predictive model
sample <- sample(seq(1, nrow(Factors_Overtime_Anemia)), replace=FALSE)
training <- Factors_Overtime_Anemia[sample[1:647],]
test <-Factors_Overtime_Anemia[sample[548:nrow(Factors_Overtime_Anemia)],]

lm_week0_training <- lm(TLAMS52 ~ Factor1_w0_flipped +  AGE + SEXCD  + ASIMPC01_A + SPLVL1, data = training)
predict_lm <- predict(lm_week0_training, data=training, newdata = test)

actuals_preds <- data.frame(cbind(actuals=test$TLAMS52, 
                                  predicteds=predict_lm))

accuracy <- data.frame(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max))
colnames(accuracy) <- "x"
accuracy <- subset(accuracy, x!= -Inf)
mean(accuracy$x)


# Severity Prediction

glm_week0_severity <- glm(Severity_w4 ~ Factor2_w0 + AGE + SEXCD + SPLVL1 + acute_duration + polytrauma + inj_update, data = Factors_Overtime_AIS, family = quasibinomial)
summary(glm_week0_severity)

library(stargazer)
stargazer(cbind(p=summary(lm_week0_Anemia_variables)$coefficients[,4],B=coef(lm_week0_Anemia_variables), confint(lm_week0_Anemia_variables)), type = "text")

lm_week0_Anemia_A <- lm(TLAMS52 ~ Factor1_w0_flipped +  AGE + SEXCD  + SPLVL1, ASIMPC01_A == "A", data = Factors_Overtime_Anemia)
lm_week0_Anemia_B <- lm(TLAMS52 ~ Factor1_w0_flipped +  AGE + SEXCD  + SPLVL1, ASIMPC01_A == "B", data = Factors_Overtime_Anemia)
lm_week0_Anemia_C <- lm(TLAMS52 ~ Factor1_w0_flipped +  AGE + SEXCD  + SPLVL1, ASIMPC01_A == "C", data = Factors_Overtime_Anemia)
lm_week0_Anemia_D <- lm(TLAMS52 ~ Factor1_w0_flipped +  AGE + SEXCD  + SPLVL1, ASIMPC01_A == "D", data = Factors_Overtime_Anemia)

summary(lm_week0_Anemia_A)
summary(lm_week0_Anemia_B)
summary(lm_week0_Anemia_C)
summary(lm_week0_Anemia_D)

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

#Another good outliers test 
library(olsrr)
ols_plot_resid_stud(lm_week0_Anemia)

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

# ncvTest and bptest are different because bptest is studendized and ncvTest is not
# => which means bptest is more robust than the original one (ncvTest)
# => go with bptest -> null is: homoscedasicity 
library(olsrr)
ols_test_breusch_pagan(lm_week0_Anemia, rhs = TRUE, multiple = TRUE)

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

#Use this to treat for heteroscadasicity or generalized additive measures or rlm?
library(lmtest)
library(sandwich)
library(caret)
cov.lm <- sqrt(diag(vcovHC(lm_week0_Anemia, type="HC")))
naive.std <- summary(lm_week0_Anemia)$coefficients[,2]
coef_table <- cbind("Estimate"=coef(lm_week0_Anemia),
                    "Naive SE" = naive.std,
                    "Robust SE" = cov.lm,
                    "p value" = summary(lm_week0_Anemia)$coefficients[,4])
library(stargazer)
stargazer(coef_table, type ="text")

#Faster code!
stargazer(cbind("Estimate"=coef(lm_week0_Anemia),
                confint(lm_week0_Anemia),
                "Naive SE" = summary(lm_week0_Anemia)$coefficients[,2],
                "Robust SE" = sqrt(diag(vcovHC(lm_week0_Anemia, type="HC"))),
                "P value" = summary(lm_week0_Anemia)$coefficients[,4]), type = "text")

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

Factors_Overtime_Anemia_noMS <- subset(Factors_Overtime_Anemia, !is.na(TLAMS52))
Factors_Overtime_Anemia_noAIS <- subset(Factors_Overtime_Anemia, !is.na(ASIMPC01_A))


UCP_RBC <- ctree(TLAMS52 ~ RBC00 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia_noMS)
plot(UCP_RBC)

UCP_HGB <- ctree(TLAMS52 ~ HGB00 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
plot(UCP_HGB)

UCP_HCT <- ctree(TLAMS52 ~ HCT00 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime_Anemia)
plot(UCP_HCT)

UCP_AIS <- ctree(ASIMPC01_A ~ RBC00*HGB00*HCT00 + SPLVL1 + AGE, data = Factors_Overtime_Anemia_noAIS)
plot(UCP_AIS)

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

Factors_Walk <- do.call("cbind", list(MS52_raw, MS00_raw, p_week0, p_week1, p_week2, p_week4, p_week8, p_week52, Bio, Other_Variables))

Factors_Walk$SEXCD <- as.factor(Factors_Walk$SEXCD) 
Factors_Walk$AGE <- as.numeric(Factors_Walk$AGE)
Factors_Walk$ASIMPC01_A <- as.factor(Factors_Walk$ASIMPC01_A)

glm_walk <- glm(Walk ~ Factor1_w0 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Walk, family = "quasibinomial")
summary(glm_walk)

drop1(glm_walk, test = "Chi")
stepAIC(glm_walk)

glm_walk_reduce <- glm(Walk ~ TLAMS00 + AGE + ASIMPC01_A + SEXCD + SPLVL1, data = Factors_Walk, family = "quasibinomial")
summary(glm_walk_reduce)

anova(glm_walk, glm_walk_reduce, test = "Chisq")

#Predictive Model 
#Split data in half to build predictive model
sample_walk <- sample(seq(1, nrow(Factors_Walk)), replace=FALSE)
training_walk <- Factors_Walk[sample_walk[1:647],]
test_walk <-Factors_Walk[sample[548:nrow(Factors_Walk)],]


glm_walk_training <- glm(Walk ~ Factor1_w0 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = training_walk, family="binomial")
predict_glm_walk <- predict(glm_walk_training, data=training_walk, newdata = test_walk, type="response")
predict_glm_walk <- as.numeric(predict_glm_walk)

walk_prediction <- ifelse(predict_glm_walk<.5, 0,1)

#Table showing predicted vs actuals
addmargins(table(walk_prediction, test_walk$Walk, dnn=c("actual", "predicted")))

#Accuracy by dividing true positive + true negative over sum
#This may overestimates because there is a dominant class -> calculate precision vs recall methods
#Precision looks at the accuracy of the positive prediction. 
#Recall is the ratio of positive instances that are correctly detected by the classifier

walk_table <-table(walk_prediction, test_walk$Walk, dnn=c("actual", "predicted"))
sum(diag(walk_table)) / sum(walk_table) #91% <- maybe overestimating

precision <- function(matrix) {
  # True positive
  tp <- matrix[2, 2]
  # false positive
  fp <- matrix[1, 2]
  return (tp / (tp + fp))
}

recall <- function(matrix) {
  # true positive
  tp <- matrix[2, 2]# false positive
  fn <- matrix[2, 1]
  return (tp / (tp + fn))
}

prec_walk <- precision(walk_table)
rec_walk <- recall(walk_table)

2 * ((prec_walk * rec_walk) / (prec_walk + rec_walk)) #=>83%


test_walk$Ordered_Walk <- factor(test_walk$Walk, levels=c("0", "1"), ordered= TRUE)
ROCpred <- roc(walk_val$Walk, walk_val$num_walk_prediction)
plot(ROCpred)

#ROC AUC vs PR AUC. if class imbalance, pick PR. if not, ROC.
#Class imbalance: yes = More of one class of data than other class
walk_val <- cbind(walk_prediction, test_walk)
walk_val <- subset(walk_val, select=c("walk_prediction", "Walk"))
walk_val$walk_prediction <- as.factor(walk_val$walk_prediction)
walk_val$num_walk_prediction <- as.numeric(walk_val$walk_prediction)

#URP for walk
Factors_Walk_noNa <- subset(Factors_Walk, !is.na(Walk))
UCP_walk <- ctree(Walk ~ Factor1_w0 + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Walk_noNa)
plot(UCP_walk)

sjp.glm(glm_walk_reduce)

#ROC for walk 
Factors_Walk$Ordered_ASIMPC01_A <- factor(Factors_Walk$ASIMPC01_A, levels=c("A", "B", "C", "D"), ordered=TRUE)
Factors_Walk$Ordered_SPLVL1 <- factor(Factors_Walk$SPLVL1, levels=c("T", "C"), ordered= TRUE)
Factors_Walk$Ordered_SEXCD <- factor(Factors_Walk$SEXCD, levels=c("1", "2"), ordered= TRUE)
Factors_Walk$Ordered_Factor1_w0 <- factor(Factors_Walk$Factor1_w0, ordered = TRUE)
Factors_Walk$Ordered_Factor2_w0 <- factor(Factors_Walk$Factor2_w0, ordered = TRUE)

W1 = roc(Walk ~ Ordered_Factor1_w0, data=Factors_Walk, smooth=TRUE)  
plot(W1, col="red", main = "ROC for Factor 1 (Blood)")

roc_wA = roc(Walk ~ Ordered_Factor1_w0, data=subset(Factors_Walk, ASIMPC01_A=="A"))  
roc_wB = roc(Walk ~ Ordered_Factor1_w0, data=subset(Factors_Walk, ASIMPC01_A=="B"))  
roc_wC = roc(Walk ~ Ordered_Factor1_w0, data=subset(Factors_Walk, ASIMPC01_A=="C"))  
roc_wD = roc(Walk ~ Ordered_Factor1_w0, data=subset(Factors_Walk, ASIMPC01_A=="D"))  

plot(roc_wA, main = "ROC of Factor 1 in Walking across each AIS grades", col="snow4", type="l", pch = 19, lty = 5)
lines(roc_wB, type="l", col = "burlywood4", pch = 16, lty = 2)
lines(roc_wC, type = "l", col = "orangered2", pch = 17, lty = 1)
lines(roc_wD, col = "tan3", type ="l", pch = 18, lty = 3)
legend(0.37,0.5, legend=c("AIS A", "AIS B", "AIS C", "AIS D"),
       lty = c(5,2,1,3),col = c("snow4", "burlywood4", "orangered2", "tan3"), 
       title = "Line Types", cex =0.8, box.lty = 0)


#Test for dispersion in glm model
library(msme)
P__disp(glm_walk)


# Does Factor 1 predict AIS grades? multinomial logistic
library(nnet)
multiple_logistic <- multinom(ASIMPC01_A ~ Factor1_w0_flipped + polytrauma + acute_duration, data = Factors_Overtime_Anemia_surg)
summary(multiple_logistic)

z <- summary(multiple_logistic)$coefficients/summary(multiple_logistic)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2

exp(coef(multiple_logistic))


AIS_grades <- polr(ASIMPC01_A ~ polytrauma + acute_duration + Factor1_w0_flipped + SEXCD + AGE + SPLVL1, data = Factors_Overtime_Anemia_surg)
AIS_table <- coef(summary(AIS_grades))
AIS_table <- cbind(AIS_table, "p value" = pnorm(abs(AIS_table[, "t value"]), lower.tail = FALSE) * 2)
exp(AIS_grades$zeta)

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


Factors_MR <- do.call("cbind", list(MS01_raw, MS00_raw, MS52_raw, p_week0, p_week1, p_week2, p_week4, p_week8, p_week52, Other_Variables, cfa_week0_only))
Factors_MR_scaled <- do.call("cbind", list(MS01_raw, MS00_raw, MS52_raw, p_week0, p_week1, p_week2, p_week4, p_week8, p_week52, Other_Variables, cfa_week0_only_scaled))

# GLM for marked Recovery (Factor 1 and 2 alone are significant but when consider other variables, only Factor 1 is sig)
Factors_MR <- Factors_MR %>% mutate_if(is.character,as.factor)
Factors_MR$SEXCD <- as.factor(Factors_MR$SEXCD) 
Factors_MR$AGE <- as.numeric(Factors_MR$AGE)
Factors_MR$Marked_Recovery <- as.factor(Factors_MR$Marked_Recovery)
Factors_MR$Factor1_w0_flipped <- Factors_MR$Factor1_w0+1
Factors_MR_noNa <- subset(Factors_MR, !is.na(Marked_Recovery))
Factors_MR_A_noNa <- subset(Factors_MR_noNa, ASIMPC01_A=="A")
Factors_MR_surg <- merge(surg_vars, Factors_MR, by.x=c("ptid"), by.y=c("PTID"))


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

glm_MR_variables <- glm(Marked_Recovery ~ RBC00*HCT00*HGB00 + AGE + SEXCD  + ASIMPC01_A + SPLVL1, data = Factors_MR_scaled, family = quasibinomial)
summary(glm_MR_variables)

drop1(glm_MR_inj, test = "Chi")

glm_MR_reduced <- glm(Marked_Recovery ~ Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = Factors_MR, family = quasibinomial)
summary(glm_MR_reduced)

glm_MR_A <- glm(Marked_Recovery ~ Factor1_w0_flipped +  AGE + SEXCD  + SPLVL1,  subset = ASIMPC01_A == "A", data = Factors_MR, family = binomial)
glm_MR_B <- glm(Marked_Recovery ~ Factor1_w0_flipped +  AGE + SEXCD  + SPLVL1,  subset = ASIMPC01_A == "B", data = Factors_MR, family = binomial)
glm_MR_C <- glm(Marked_Recovery ~ Factor1_w0_flipped +  AGE + SEXCD  + SPLVL1,  subset = ASIMPC01_A == "C", data = Factors_MR, family = binomial)
glm_MR_D <- glm(Marked_Recovery ~ Factor1_w0_flipped +  AGE + SEXCD  + SPLVL1,  subset = ASIMPC01_A == "D", data = Factors_MR, family = binomial)

summary(glm_MR_A)
summary(glm_MR_B)
summary(glm_MR_C)
summary(glm_MR_D)

anova(glm_MR_reduced, glm_MR, test = "Chisq")

#Prediction values of glm: 
glm_MR_probs <- predict(glm_MR, type="response")
summary(glm_MR_probs) #pick the mean?
glm_MR_pred <- ifelse(glm_MR_probs>0.13, 1, 0) 
table(glm_MR_pred, Factors_MR_noNa$Marked_Recovery)
mean(glm_MR_pred == Factors_MR_noNa$Marked_Recovery)

#Prediction values of glm based on AIS-A alone: 
glm_MR_A_probs <- predict(glm_MR_A, type="response")
glm_MR_A_pred <- ifelse(glm_MR_A_probs>0.14, 1, 0)
table(glm_MR_A_pred, Factors_MR_A_noNa$Marked_Recovery)
mean(glm_MR_A_pred == Factors_MR_A_noNa$Marked_Recovery)

#test for dispersion 
#One way to see if model is overdispersed is looking at Residual Deviance : degrees of freedom
# more than 1 -> Overdispersion
library(msme)
P__disp(glm_MR)

#Ctree

UCP_MR <- ctree(Marked_Recovery ~ Factor1_w0 + SPLVL1 + AGE + ASIMPC01_A, data = subset(Factors_MR, !is.na(Marked_Recovery)))
plot(UCP_MR)

#Model diagnostics - Test for linearity and multicollinearity 
crPlots(glm_MR)
vif(glm_MR)

library(boot)
glm.diag.plots (glm_MR)

# ROC for Marked Recovery 
library(pROC) #Only works with numeric or ordered factors so if it is factors, have to order it
Factors_MR$Ordered_ASIMPC01_A <- factor(Factors_MR$ASIMPC01_A, levels=c("A", "B", "C", "D"), ordered=TRUE)
Factors_MR$Ordered_SPLVL1 <- factor(Factors_MR$SPLVL1, levels=c("T", "C"), ordered= TRUE)
Factors_MR$Ordered_SEXCD <- factor(Factors_MR$SEXCD, levels=c("1", "2"), ordered= TRUE)
Factors_MR$Ordered_Factor1_w0 <- factor(Factors_MR$Factor1_w0_flipped, ordered = TRUE)
Factors_MR$Ordered_Factor2_w0 <- factor(Factors_MR$Factor2_w0, ordered = TRUE)

AUC = roc(Marked_Recovery ~ Factor1_w0 + Ordered_ASIMPC01_A + Ordered_SPLVL1 + AGE + Ordered_SEXCD, data=Factors_MR, plot = TRUE, percent = TRUE, ci = TRUE) 

F1 = roc(Marked_Recovery ~ Ordered_Factor1_w0, data=Factors_MR)  
plot(F1, col="red", main = "ROC for Factor 1 (Blood)")

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

roc_A = roc(Marked_Recovery ~ Ordered_Factor1_w0, data=subset(Factors_MR, ASIMPC01_A=="A"), smooth=TRUE)  
roc_B = roc(Marked_Recovery ~ Ordered_Factor1_w0, data=subset(Factors_MR, ASIMPC01_A=="B"), smooth=TRUE)  
roc_C = roc(Marked_Recovery ~ Ordered_Factor1_w0, data=subset(Factors_MR, ASIMPC01_A=="C"), smooth=TRUE)  
roc_D = roc(Marked_Recovery ~ Ordered_Factor1_w0, data=subset(Factors_MR, ASIMPC01_A=="D"))  

plot(roc_A, main = "ROC of Factor 1 in each AIS grades", col="snow4", type="l", pch = 19, lty = 5)
lines(roc_B, type="l", col = "burlywood4", pch = 16, lty = 2)
lines(roc_C, type = "l", col = "orangered2", pch = 17, lty = 1)
lines(roc_D, col = "tan3", type ="l", pch = 18, lty = 3)
legend(0.37,0.5, legend=c("AIS A", "AIS B", "AIS C", "AIS D"),
       lty = c(5,2,1,3),col = c("snow4", "burlywood4", "orangered2", "tan3"), 
      title = "Line Types", cex =0.8, box.lty = 0)


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
cbind(Odds=exp(coef(glm_MR)), exp(confint(glm_MR)),p=summary(glm_MR)$coefficients[,4])

Variables_glm <- cbind(Odds=exp(coef(glm_MR_variables)), exp(confint(glm_MR_variables)), p=summary(glm_MR_variables)$coefficients[,4])
stargazer(Variables_glm, type = "text")

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
Factors_Overtime_RBC_drugs$Factor1_w0_flipped <- Factors_Overtime_RBC_drugs$Factor1_w0+1

lm_RBC_drugs <- lm(TLAMS52 ~ RBC_Drug, data = Factors_Overtime_RBC_drugs)
lm_factor1 <- lm(TLAMS52 ~ Factor1_w0_flipped + AGE + SPLVL1 + ASIMPC01_A + RBC_Drug + SEXCD, data = Factors_Overtime_RBC_drugs)
lm_reduced_factor1 <- lm(TLAMS52 ~ Factor1_w0_flipped + AGE + SPLVL1 + ASIMPC01_A + RBC_Drug, data = Factors_Overtime_RBC_drugs)
lm_RBC_factor1 <- lm(TLAMS52 ~ Factor1_w0_flipped + RBC_Drug, data = Factors_Overtime_RBC_drugs)
anova(lm_RBC_factor1, lm_factor1, test = "Chisq")

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

#Factor 1 vs MS
ggplot(data = subset(Factors_Overtime_Anemia, !is.na(Severity)),aes(x=TLAMS52, y=Factor1_w0_flipped))+
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  ggtitle("Blood factor across severities in LEMS 52") +
  facet_grid(. ~ Severity, scales = "free_x")


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
Factors_Overtime_Death_Sur <- do.call("cbind", list(Factors_Overtime_Death_Sur, cfa_week0_only)) 

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

Factors_Overtime_Death_Sur$Age_binary <- ifelse(Factors_Overtime_Death_Sur$AGE <= 50, 0, 1)
Factors_Overtime_Death_Sur$Age_binary <- as.factor(Factors_Overtime_Death_Sur$Age_binary)  

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
Hazard_Factor3 <- coxph(Surv(Survival_days, Survival_Rates==0)~ Factor3_w0_flipped + AGE + SEXCD + ASIMPC01_A + SPLVL1, data = Factors_Overtime_Death_Sur)
summary(Hazard_Factor3)

Hazard_Factor4 <- coxph(Surv(Survival_days, Survival_Rates==0)~ Factor4_w0_flipped + AGE + SEXCD + ASIMPC01_A + SPLVL1, data = Factors_Overtime_Death_Sur)
summary(Hazard_Factor4)

Hazard_Factors <- coxph(Surv(Survival_days, Survival_Rates==0)~ Factor3_w0_flipped + Factor4_w0_flipped + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Death_Sur)
summary(Hazard_Factors)

Hazard_Variables_Factor3 <- coxph(Surv(Survival_days, Survival_Rates==0)~ CK000*ALT00*AST00 + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Death_Sur)
summary(Hazard_Variables_Factor3)

Hazard_Variables_Factor4 <- coxph(Surv(Survival_days, Survival_Rates==0)~ BUN00*BC900*BUA00*P0400 + AGE + SEXCD + SPLVL1 + ASIMPC01_A, data = Factors_Overtime_Death_Sur)
summary(Hazard_Variables_Factor4)

Hazard_Factor3_age <- coxph(Surv(Survival_days, Survival_Rates==0)~ Factor3_w0_flipped + Age_binary, data = Factors_Overtime_Death_Sur)
summary(Hazard_Factor3_age)

#37 out of 47 dead patients have AIS-A
Hazard_Factor3_A <- coxph(Surv(Survival_days, Survival_Rates==0)~ Factor3_w0_flipped + AGE + SEXCD + SPLVL1, subset = ASIMPC01_A == "A", data = Factors_Overtime_Death_Sur)
Hazard_Factor4_A <- coxph(Surv(Survival_days, Survival_Rates==0)~ Factor4_w0_flipped + AGE + SEXCD + SPLVL1, subset = ASIMPC01_A == "A", data = Factors_Overtime_Death_Sur)

summary(Hazard_Factor3_A)
summary(Hazard_Factor4_A)

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
Wei_Factor3 <- survreg(Surv(Survival_days, Survival_Rates == 0) ~ Factor3_w0_flipped + AGE + SEXCD + SPLVL1 + ASIMPC01_A, dist = "weibull", data = Factors_Overtime_Death_Sur)
Converted_Wei <- ConvertWeibull(Wei_Factor3, conf.level = 0.95)

Wei_Factor4 <- survreg(Surv(Survival_days, Survival_Rates == 0) ~ Factor4_w0_flipped + AGE + SEXCD + SPLVL1 + ASIMPC01_A, dist = "weibull", data = Factors_Overtime_Death_Sur)
Converted_Wei4 <- ConvertWeibull(Wei_Factor4, conf.level = 0.95)

Wei_Variables_Factor3 <- survreg(Surv(Survival_days, Survival_Rates == 0) ~ CK000*ALT00*AST00 + AGE + SEXCD + SPLVL1 + ASIMPC01_A, dist = "weibull", data = Factors_Overtime_Death_Sur)
Converted_Wei_Variables3 <- ConvertWeibull(Wei_Variables_Factor3, conf.level = 0.95)

Wei_Variables_Factor4 <- survreg(Surv(Survival_days, Survival_Rates == 0) ~ BUN00*BC900*BUA00*P0400 + AGE + SEXCD + SPLVL1 + ASIMPC01_A, dist = "weibull", data = Factors_Overtime_Death_Sur)
Converted_Wei_Variables4 <- ConvertWeibull(Wei_Variables_Factor4, conf.level = 0.95)

Wei_Factor3_A <- survreg(Surv(Survival_days, Survival_Rates==0)~ Factor3_w0_flipped + AGE + SEXCD + SPLVL1, subset = ASIMPC01_A == "A",dist = "weibull", data = Factors_Overtime_Death_Sur)
Converted_Wei_Factor3_A <- ConvertWeibull(Wei_Factor3_A, conf.level = 0.95)

Wei_Factor4_A <- survreg(Surv(Survival_days, Survival_Rates==0)~ Factor4_w0_flipped + AGE + SEXCD + SPLVL1, subset = ASIMPC01_A == "A",dist = "weibull", data = Factors_Overtime_Death_Sur)
Converted_Wei_Factor4_A <- ConvertWeibull(Wei_Factor4_A, conf.level = 0.95)


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
forest_model(Hazard_Factor3)
forest_model(Hazard_Factor4)

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

Cor_cfa_52 <- cor(cfa_week52, use="pairwise.complete.obs")
Cor_cfa_52 <- as.matrix(Cor_cfa_52)
corrplot(Cor_cfa_52, method = "color", type = "upper", order = "hclust",
         addCoef.col = "black", number.cex=0.35)

RBC_cor <- as.matrix(cor(RBC, use="pairwise.complete.obs"))
colnames(RBC_cor) <- c("Red Blood Cells", "Hematocrit", "Hemoglobin")
rownames(RBC_cor) <- c("Red Blood Cells", "Hematocrit", "Hemoglobin")
corrplot(RBC_cor, method = "color", type = "upper", order = "hclust")

RBC_corplot <- ggplot(Factors_Overtime_Anemia_scaled, aes(x=TLAMS52, y=RBC00))+
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("text", x=15, y=5, label="r = 0.2, p < 0.05", colour = "red") +
  theme_bw() +
  xlab("LEMS at Week 52") +
  ylab("Red Blood Cells Count") + 
  ggtitle("Red Blood Cells Count")

ggplot(data=subset(Factors_Overtime_Anemia, !is.na(ASIMPC01_A)),aes(x=TLAMS52, y=RBC00))+
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("LEMS at Week 52") +
  ylab("Red Blood Cells Count") + 
  ggtitle("Red Blood Cells Count Stratified") +
  facet_grid(. ~ ASIMPC01_A, scale = "free_x")

HCT_corplot <- ggplot(Factors_Overtime_Anemia, aes(x=TLAMS52, y=HCT00))+
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("text", x=18, y=5, label="r = 0.21, p < 0.05", colour = "red") +
  theme_bw() +
  xlab("LEMS at Week 52") +
  ylab("Hematocrit") + 
  ggtitle("Hematocrit")

ggplot(data=subset(Factors_Overtime_Anemia_scaled, !is.na(ASIMPC01_A)),aes(x=TLAMS52, y=HCT00))+
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("LEMS at Week 52") +
  ylab("Hematocrit") + 
  ggtitle("Hematocrit Stratified") +
  facet_grid(. ~ ASIMPC01_A, scale = "free_x")

HGB_corplot <- ggplot(Factors_Overtime_Anemia_scaled, aes(x=TLAMS52, y=HGB00))+
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("text", x=18, y=5, label="r = 0.21, p < 0.05", colour = "red") +
  theme_bw() +
  xlab("LEMS at Week 52") +
  ylab("Hemoglobin") + 
  ggtitle("Hemoglobin")

ggplot(data=subset(Factors_Overtime_Anemia_scaled, !is.na(ASIMPC01_A)),aes(x=TLAMS52, y=HGB00))+
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("LEMS at Week 52") +
  ylab("Hemoglobin") + 
  ggtitle("Hemoglobin Stratified") +
  facet_grid(. ~ ASIMPC01_A, scale = "free_x")

Factor1_corplot <- ggplot(Factors_Overtime_Anemia, aes(x=TLAMS52, y=Factor1_w0_flipped))+
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("text", x=18, y=5, label="r = 0.22, p < 0.05", colour = "red") +
  theme_bw() +
  xlab("LEMS at Week 52") +
  ylab("Factor 1 (Blood)") + 
  ggtitle("Factor 1 (Blood)")

multiplot(RBC_corplot, HCT_corplot, HGB_corplot, Factor1_corplot, cols = 2)

#Clustering of AIS grades
#Using k-means methods
clust = Factors_Overtime_Anemia[,c("ASIMPC01_A","HGB00","HCT00","RBC00")]
clust<-clust[complete.cases(clust),]

clus2 <- clust
clus2$ASIMPC01_A <- NULL
cluster = kmeans(clus2, 4)

plot(clus2, col = (cluster$cluster+1), pch=20, cex=2)

table(clust$ASIMPC01_A, cluster$cluster)

#using k-means methods
pam.result <- pam(clus2, 4)
table(pam.result$clustering, clust$ASIMPC01_A)
plot(pam.result)

plot(clust[c("HGB00", "TLAMS52")], col = cluster$cluster)
points(cluster$centers[, c("HGB00", "TLAMS52")],
       col = 1:3, pch = 8, cex = 2)

ggplot(data=subset(Factors_Overtime_Anemia, !is.na(ASIMPC01_A)), aes(x=TLAMS52, y=HGB00, colour=ASIMPC01_A))+
  geom_boxplot()+
  geom_jitter(width=0.2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  

describeBy(Factors_Overtime_Anemia$HGB00, Factors_Overtime_Anemia$ASIMPC01_A)

#try multiple factor analysis
library(FactoMineR)
library(factoextra)

mfa <- MFA(clust, group = c(1,3),
           type = c("n", "s"),
           name.group = c("AIS", "Blood"),
           graph = FALSE)
fviz_mfa_var(mfa, "group")

fviz_ellipses(mfa, c("ASIMPC01_A"), repel = TRUE)

#try linear discriminant analysis 
library(MASS)
lda <- lda(ASIMPC01_A ~ RBC00+HCT00+HGB00, data=clust)
plot(lda, col = as.integer(clust$ASIMPC01_A))

lda.pred <- predict(lda)
clust$lda <- lda.pred$class

table(clust$lda, clust$ASIMPC01_A)

library(klaR)
partimat(ASIMPC01_A ~ HCT00+HGB00, data=clust, method="qda")

#Changing to time format 
surgical$ctsp_moment <- as.POSIXct(surgical$ctsp_moment, format='%Y/%m/%d %H:%M')
surgical$ctst_moment <- as.POSIXct(surgical$ctst_moment, format='%Y/%m/%d %H:%M')
surgical$otsp_moment <- as.POSIXct(surgical$otsp_moment, format='%Y/%m/%d %H:%M')
surgical$otst_moment <- as.POSIXct(surgical$otst_moment, format='%Y/%m/%d %H:%M')
surgical$pasp_moment <- as.POSIXct(surgical$pasp_moment, format='%Y/%m/%d %H:%M')
surgical$past_moment <- as.POSIXct(surgical$past_moment, format='%Y/%m/%d %H:%M')
surgical$surgdt_update <- as.POSIXct(surgical$surgdt_update, format='%Y/%m/%d')
surgical$emc_moment <- as.POSIXct(surgical$emc_moment, format='%Y/%m/%d %H:%M')
surgical$aracut_moment <- as.POSIXct(surgical$aracut_moment, format='%Y/%m/%d %H:%M')
surgical$diacut_moment <- as.POSIXct(surgical$diacut_moment, format='%Y/%m/%d %H:%M')
surgical$mpsssp_moment <- as.POSIXct(surgical$mpsssp_moment, format='%Y/%m/%d %H:%M')
surgical$mpssst_moment <- as.POSIXct(surgical$mpssst_moment, format='%Y/%m/%d %H:%M')

ID_MS <- c("TLAMS52", "PTID")
ID_MS <- Factors_Overtime_Anemia[ID_MS]

surgical$cts_duration <- difftime(surgical$ctsp_moment, surgical$ctst_moment,units="hours")
surgical$cts_duration <- as.numeric(surgical$cts_duration)
surgical$cts_duration[surgical$cts_duration<0] = NA

surgical$ots_duration <- difftime(surgical$otsp_moment, surgical$otst_moment, units="hours")
surgical$ots_duration <- as.numeric(surgical$ots_duration)
surgical$ots_duration[surgical$ots_duration<0] = NA

surgical$pas_duration <- difftime(surgical$pasp_moment, surgical$past_moment, units="hours")
surgical$pas_duration <- as.numeric(surgical$pas_duration)
surgical$pas_duration[surgical$pas_duration<0] = NA

surgical$emc_acute_duration <- difftime(surgical$aracut_moment, surgical$emc_moment, units="hours")
surgical$emc_acute_duration <- as.numeric(surgical$emc_acute_duration)
surgical$emc_acute_duration[surgical$emc_acute_duration<0] = NA

surgical$acute_duration <- difftime(surgical$diacut_moment, surgical$aracut_moment, units="hours")
surgical$acute_duration <- round(as.numeric(surgical$acute_duration), digits=0)
surgical$acute_duration[surgical$acute_duration<0] = NA

surgical$mpss_duration <- difftime(surgical$mpsssp_moment, surgical$mpssst_moment, units="hours")
surgical$mpss_duration <- as.numeric(surgical$mpss_duration)
surgical$mpss_duration[surgical$mpss_duration<0] = NA

surgical$wait_duration <- difftime(surgical$aracut_moment, surgical$emc_moment, units="hours")
surgical$wait_duration <- as.numeric(surgical$wait_duration)
surgical$wait_duration[surgical$wait_duration<0] = NA

surgical_MS <- merge(surgical, ID_MS, by.x=c("ptid"), by.y=c("PTID"))
surgical_MS <- merge(surgical_MS, etio, by.x=c("ptid"), by.y=c("ptid"))

surgical_MS_noNa <- surgical_MS[complete.cases(surgical_MS[ , 118:122, 85]),]

surgical_MS$inj_comb <- rowSums(surgical_MS[,c("cardcd01", "cpulcd01", "eentcd01", 
                                               "gastcd01", "genicd01", "headcd01" , 
                                               "musccd01", "pulmcd01", "skincd01")])

surgical_MS$polytrauma <- ifelse(surgical_MS$inj_comb < 3, "0", "1")

surg_vars <- c("ptid", "acute_duration", "polytrauma", "inj_update")
surg_vars <- surgical_MS[surg_vars]

surg_other <- lm(TLAMS52 ~ bpsys3 + bpdia3 + heart3 + resprt3 + tempc3 + acute_duration + wait_duration + polytrauma, data=surgical_MS)
backward(surg_other)

surg <- lm(TLAMS52 ~ bpsys3 + heart3 + acute_duration, data = surgical_MS_noNa)

anova(surg, surg_other)

#Transform into long format in Factors data to see the trajectory of Factor 1
library(reshape)
Factors_Overtime_long <- melt(Factors_MR, id.vars=c("TLAMS52", "PTID", "SEXCD", "AGE", "ASIMPC01_A", "Marked_Recovery", "SPLVL1"),
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
Factors_Overtime_long <- separate(data = Factors_Overtime_long, col = variable, 
                             into = c("Factors", "TimeAfterInjury"), sep = "\\_")

Factors_Overtime_long$TimeAfterInjury = factor(Factors_Overtime_long$TimeAfterInjury, levels=c('w0','w1','w2','w4','w8','w52'))
Factors_Overtime_long$TimeAfterInjury_numeric <- gsub("w", "", Factors_Overtime_long$TimeAfterInjury)
Factors_Overtime_long$TimeAfterInjury_numeric <- as.numeric(Factors_Overtime_long$TimeAfterInjury_numeric)
Factors_Overtime_long$TimeAfterInjury_numeric_scaled <- scale(Factors_Overtime_long$TimeAfterInjury_numeric)

Factors_Overtime_long_noW2 <- subset(Factors_Overtime_long, TimeAfterInjury_numeric!="8" & TimeAfterInjury_numeric!="52")
#Factor 1 changes over time for each patient
xyplot(value~TimeAfterInjury|PTID, col.line = "red", layout = c(9,5), ylim=c(-2,2),
       grid = TRUE, type = c("p", "r"), data = subset(Factors_Overtime_long_noW2, Factors=="Factor1"))

#For the error term, anything that is to the right of "|" denotes grouping variable (error term)
# anything that is to the left of "|" denotes term varies between groups (1:random intercept, +X: random slope)
Factor_time <- lmer(value ~ TimeAfterInjury_numeric+ (1+TimeAfterInjury_numeric|PTID), REML = FALSE,
            data = subset(Factors_Overtime_long_noW2,
                          Factors=="Factor2"))

summary(Factor_time)
coefs_factor <- coef(Factor_time)$PTID[1:2]
colnames(coefs_factor) <- c("intercept", "slope")
library(data.table)
setDT(coefs_factor, keep.rownames = TRUE)[]
colnames(coefs_factor)[colnames(coefs_factor)=="rn"] <- "PTID"

Factors_Overtime_coefs <- merge(coefs_both, Factors_MR, all=TRUE)
Factors_Overtime_coefs_update <- merge(Factors_Overtime_coefs, MS_only, by="PTID")
Factors_Overtime_coefs_update$Marked_Recovery <- as.factor(Factors_Overtime_coefs_update$Marked_Recovery)


ggplot(data=coefs_factor, aes(x=slope, y=intercept))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE)

#Factor 1 alone is not significant over time
ggplot(data=subset(Factors_Overtime_long_noW2, !is.na(ASIMPC01_A)), aes(x=TimeAfterInjury, y=value))+
  geom_point(aes(group=PTID))+
  geom_line(aes(group=PTID))+
  facet_grid(.~ASIMPC01_A)

plot(Effect(c("TimeAfterInjury_numeric","ASIMPC01_A"),Factor_time))

plot(Factor_time)
abline(h=0,lwd=2)

qqnorm(resid(Factor_time))
qqline(resid(Factor_time))

qqnorm(ranef(Factor_time)$PTID[,1])
qqline(ranef(Factor_time)$PTID[,1])

#The first 50
Factors_50 <- head(Factors_MR, 50)
Factors_50_long <- reshape(Factors_50, varying = c("Factor1_w0", "Factor1_w1", "Factor1_w2","Factor1_w4", "Factor1_w8", "Factor1_w52"),
                                    v.names = "Factors", timevar = "TimeAfterInjury", times = c("w0", "w1","w2", "w4", "w8", "w52"), direction = "long")

Factors_50_long$TimeAfterInjury = factor(Factors_50_long$TimeAfterInjury, levels=c('w0','w1','w2','w4','w8','w52'))
Factors_50_long$TimeAfterInjury_numeric <- gsub("w", "", Factors_50_long$TimeAfterInjury)
Factors_50_long$TimeAfterInjury_numeric <- as.numeric(Factors_50_long$TimeAfterInjury_numeric)
Factors_50_long$TimeAfterInjury_numeric_scaled <- scale(Factors_50_long$TimeAfterInjury_numeric)


plot(ctree(Marked_Recovery ~ SPLVL1 + SEXCD + AGE + ASIMPC01_A + Factors, data = subset(Factors_Overtime_coefs_update, !is.na(Severity))))



ggplot(data = subset(Factors_50_long, TimeAfterInjury!="w8"&TimeAfterInjury!="w52"), aes(x=TimeAfterInjury, y=Factors, group=PTID)) +
  geom_point()+
  geom_line()+
  geom_smooth(method="lm",formula = y~poly(x,3, raw=TRUE), se=FALSE,colour="red")
  
  facet_grid(.~ASIMPC01_A, space="free")
 

#Interacted with AIS grades, they are significant 
ggplot(data=subset(Factors_Overtime_long_noW2, !is.na(ASIMPC01_A)), aes(x=TimeAfterInjury, y=value))+
  geom_boxplot(aes(fill = ASIMPC01_A), width = 2)+
  coord_cartesian(ylim=c(-2, 2)) +
  scale_fill_manual(values = c("lightgoldenrod1", "darkgoldenrod", "paleturquoise", "paleturquoise4")) +
  ggtitle("Factor 1 Values Over Time")+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 12))+
  facet_grid(.~TimeAfterInjury, scales="free_x")

#Transform into long format for Motor scores trajectory  
Factors_MS_Overtime <- do.call("cbind", list(MS00_raw, MS01_raw, MS04_raw, MS08_raw, MS16_raw, MS26_raw, MS52_raw, 
                                             p_week0, p_week1, p_week2, p_week4, p_week8, p_week52, 
                                             Other_Variables, cfa_week0_only))

Factors_MS_Overtime_long <- reshape(Factors_MS_Overtime, varying = c("TLAMS00", "TLAMS01", "TLAMS04", "TLAMS08","TLAMS16","TLAMS26", "TLAMS52"),
                                       v.names = "LEMS", timevar = "TimeAfterInjury", times = c("w0", "w1", "w4", "w8","w16", "w26","w52"), direction = "long")

Factors_MS_Overtime_long$TimeAfterInjury = factor(Factors_MS_Overtime_long$TimeAfterInjury, levels=c('w0','w1','w4','w8','w16','w26','w52'))
Factors_MS_Overtime_long$TimeAfterInjury_numeric <- gsub("w", "", Factors_MS_Overtime_long$TimeAfterInjury)
Factors_MS_Overtime_long$TimeAfterInjury_numeric <- as.numeric(Factors_MS_Overtime_long$TimeAfterInjury_numeric)
Factors_MS_Overtime_long$TimeAfterInjury_numeric_scaled <- scale(Factors_MS_Overtime_long$TimeAfterInjury_numeric)

Factors_MS_Overtime_long_noW <- subset(Factors_MS_Overtime_long, TimeAfterInjury_numeric!="8" & TimeAfterInjury_numeric!="52" & TimeAfterInjury_numeric!="16" & TimeAfterInjury_numeric!="26")

#Each individual 
xyplot(LEMS~TimeAfterInjury|PTID, col.line = "red", layout = c(9,5),
       grid = TRUE, type = c("p", "r"), data = Factors_MS_Overtime_long)

MS_time <- lmer(LEMS ~ TimeAfterInjury_numeric+(1+TimeAfterInjury_numeric|PTID), REML=FALSE,
                    data = Factors_MS_Overtime_long_noW)
summary(MS_time)

coefs_MS <- coef(MS_time)$PTID[1:2]
colnames(coefs_MS) <- c("intercept", "slope")

setDT(coefs_MS, keep.rownames = TRUE)[]
colnames(coefs_MS)[colnames(coefs_MS)=="rn"] <- "PTID"


ggplot(data=coefs_MS, aes(x=intercept, y=slope))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE)

# MS overtime alone is not sig
ggplot(data=subset(Factors_MS_Overtime_long, !is.na(ASIMPC01_A)), aes(x=TimeAfterInjury, y=LEMS))+
  geom_boxplot()

# But it is when interacted with AIS grades
ggplot(data=subset(Factors_MS_Overtime_long, !is.na(ASIMPC01_A)), aes(x=TimeAfterInjury, y=LEMS))+
  geom_boxplot(aes(fill = ASIMPC01_A), width = 2)+
  scale_fill_manual(values = c("lightgoldenrod1", "darkgoldenrod", "paleturquoise", "paleturquoise4")) +
  ggtitle("MS Values Over Time")+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 12))+
  facet_grid(.~TimeAfterInjury, scales="free_x")

coefs_both <- merge(coefs_factor, coefs_MS, all=TRUE, by="PTID")
colnames(coefs_both) <- c("PTID", "intercept_Factor", "slope_Factor", 
                           "intercept_MS", "slope_MS")


ggplot(data=coefs_both, aes(x=intercept_Factor, y=slope_MS))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE)


ggplot(data=Factors_MS_Overtime_long, aes(x=TimeAfterInjury, y=LEMS, group=PTID))+
  geom_point()+
  geom_line()+
  geom_smooth(aes(group = 1, colour = "Trendline"),method="loess",formula = y~x, se=FALSE)

#Drugs?
Drugs <- do.call("cbind", list(Other_Variables,MS00_raw, MS01_raw, MS04_raw, 
                               MS08_raw, MS16_raw, MS26_raw, MS52_raw,
                                cfa_week0_only, cfa_week1, cfa_week2, 
                                cfa_week4, cfa_week8, cfa_week52))

Anh_1 <- as.character(Anh_1$ptid)
Drugs_1<-Drugs[Drugs$PTID %in% Anh_1,]



Drug_long <-melt(Drugs, id.vars=c("PTID", "SEXCD", "AGE", "ASIMPC01_A","SPLVL1"),
                                           measure.vars=c("HGB00", "HGB01", "HGB02",
                                                          "HGB04", "HGB08","HGB52"))

Drug_long$variable <- gsub("HGB", "w", Drug_long$variable)
Drug_long$variable = factor(Drug_long$variable, levels=c('w00','w01','w02','w04','w08','w52'))
setnames(Drug_long, c("variable","value"), c("Time", "Values"))
Drug_long$Time_numeric <- gsub("w", "", Drug_long$Time)
Drug_long$Time_numeric <- as.numeric(Drug_long$Time_numeric)
Drug_long$Time_numeric_scaled <- scale(Drug_long$Time_numeric)


Blood_lmer <-lmer(Values ~ Time_numeric_scaled+(1+Time_numeric_scaled|PTID), REML = FALSE,
                    data = Drug_long)

coefs_blood <- coef(Blood_lmer)$PTID[1:2]
colnames(coefs_blood) <- c("intercept", "slope")

setDT(coefs_blood, keep.rownames = TRUE)[]
colnames(coefs_blood)[colnames(coefs_blood)=="rn"] <- "PTID"



xyplot(value~variable|PTID, col.line = "red", layout = c(9,5),
       grid = TRUE, type = c("p", "r"), data = Drug_long)



ggplot(data = Drug_long, aes(x=Time, y=Values, group=PTID)) +
  geom_point()+
  geom_line()+
  geom_smooth(aes(group = 1, colour = "Trendline"),method="lm",formula = y~poly(x,3), se=FALSE)



