cor(cfa1_week0, use="pairwise.complete.obs")

#Scaling data of week 1
cfa1_week1 <- scale(cfa_week1, center = TRUE, scale = TRUE)
cfa1_week1 <- data.frame(cfa1_week1)
cfa1_week1

Others <- subset(Other_Variables, select=MAP:INFDOSE1)
Others <- as.data.frame(Others)


cfa_week1 <- cbind(MS52_raw, cfa_week1)
Categorical <- cbind(cfa1_week1, Other_Variables)



library(psych)
scree <- fa.parallel(cfa1_week1, fm = "ml", fa = "fa")
scree



fa <- fa(cfa1_week1[-1], nfactors=6, rotate="varimax", fm="ml")
fa.diagram(fa)

print(fa$loadings, cutoff = 0.5)

#Keep original factors of week 0 to see if factor 1 changes across time
model1 <- 'Factor 1 =~ HGB01 + HCT01 + RBC01 
Factor 2 =~ ALB01 + CAB01 + CHO01 + STP01
Factor 3 =~ CK001 + ALT01 + AST01
Factor 4 =~ BUN01 + BC901 + BUA01 + P0401
Factor 5 =~ MCV01 + MCH01
Factor 6 =~ DBN01 + DDN01 + WBC01 + DLN01 + DJN01
AST01 ~~ 0*AST01
MCV01 ~~ 0*MCV01
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

#Scaling data for week 2
cfa1_week2 <- scale(cfa_week2, center = TRUE, scale = TRUE)
cfa1_week2 <- data.frame(cfa1_week2)
cfa1_week2

cfa_week2 <- cbind(MS52_raw, cfa_week2)


library(psych)
scree <- fa.parallel(cfa1_week2[-1], fm = "ml", fa = "fa")
scree

fa <- fa(cfa1_week2[-1], nfactors=6, rotate="varimax", fm="ml")
fa.diagram(fa)

#Keep original factors of week 0 to see if factor 1 changes across time

model2 <- 'Factor 1 =~ HGB02 + HCT02 + RBC02 
           Factor 2 =~ STP02 + ALB02 + CAB02 + CHO02
           Factor 3 =~ AST02 + CK002 + ALT02 
           Factor 4 =~ BUN02 + BC902 + BUA02 + P0402
           Factor 5 =~ MCV02 + MCH02
           Factor 6 =~ DBN02 + DDN02 + WBC02 + DLN02 + DJN02'

f1_week2 <- sem(model2, data = cfa1_week2, std.lv=TRUE, na.action=na.omit, missing = "fiml")
summary(f1_week2, standardized = TRUE, fit.measures = TRUE)


semPaths(f1_week2,"std", whatLabels="std", intercepts=FALSE, style="lisrel", residuals = FALSE,
         nCharNodes=0, 
         nCharEdges=0,
         curveAdjacent = TRUE,title=TRUE, layout="tree2",curvePivot=FALSE)

p_week2 <- lavPredict(f1_week2)
p_week2 <- as.data.frame(p_week2)

colnames(p_week2) <- c("Factor1_w2", "Factor2_w2", "Factor3_w2", "Factor4_w2", "Factor5_w2")

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
           Factor 5 =~ MCV04 + MCH04'

f1_week4 <- sem(model4, data = cfa1_week4, std.lv=TRUE, na.action=na.omit, missing = "fiml")
summary(f1_week4, standardized = TRUE, fit.measures = TRUE)

semPaths(f1_week4,"std", whatLabels="std", intercepts=FALSE, style="lisrel", residuals = FALSE,
         nCharNodes=0, 
         nCharEdges=0,
         curveAdjacent = TRUE,title=TRUE, layout="tree2",curvePivot=FALSE)

p_week4 <- lavPredict(f1_week4)
p_week4 <- as.data.frame(p_week4)

colnames(p_week4) <- c("Factor1_w4", "Factor2_w4", "Factor3_w4", "Factor4_w4", "Factor5_w4")

#Raw and Scaled scores for MS52
MS52_scaled <- as.data.frame(scale(cfa_week0[1], scale = TRUE, center = TRUE))
MS52_raw <- cfa_week0[1]

#Combined dataframe
Factors_Overtime <- do.call("cbind", list(MS52_scaled, p_week0, p_week1, p_week2, p_week4, Other_Variables))

#Transforming data to long format

#Taking Correlation r values between MS52 and Factors across time and put into df
Cor_Factors <- cor(Factors_Overtime[1:21], use="pairwise.complete.obs")
Cor_Factors <- head(Cor_Factors, 1)
Cor_Factors <- as.data.frame(Cor_Factors)

#Transform into long format
library(reshape)
Cor_Factors_long <- melt(Cor_Factors, id.vars=c("TLAMS52"),
                              measure.vars=c("Factor1_w0", "Factor1_w1", "Factor1_w2", "Factor1_w4",
                                             "Factor2_w0", "Factor2_w1", "Factor2_w2", "Factor2_w4",
                                             "Factor3_w0", "Factor3_w1", "Factor3_w2", "Factor3_w4",
                                             "Factor4_w0", "Factor4_w1", "Factor4_w2", "Factor4_w4",
                                             "Factor5_w0", "Factor5_w1", "Factor5_w2", "Factor5_w4"),
                              variable.name ="Time After Injury", 
                              value.name ="Factors")

#Splitting into 2 separate columns
library(tidyr)
Cor_Factors_long <- separate(data = Cor_Factors_long, col = variable, 
                                  into = c("Factors", "TimeAfterInjury"), sep = "\\_")

#Plot across time
library(ggplot2)
p <- ggplot(Cor_Factors_long, aes(x=TimeAfterInjury, y=value)) +
  geom_jitter()+
  geom_smooth(color='red', method='lm', formula=y~x, aes(group=1), se=FALSE) +
  ggtitle("Correlation between Motor Scores at week 52 and Factors from week 0 to week 4")+
  theme_bw() +
  labs(x="Time After Injury", y="r values at MS week 52") +
  facet_grid(.~Factors)
p

#Building UCP

Factors_long <- melt(Factors_Overtime, id.vars=c("TLAMS52", "AGE", "SEXCD", "ASIMPC01_A", "SPLVL1"),
                         measure.vars=c("Factor1_w0", "Factor1_w1", "Factor1_w2", "Factor1_w4",
                                        "Factor2_w0", "Factor2_w1", "Factor2_w2", "Factor2_w4",
                                        "Factor3_w0", "Factor3_w1", "Factor3_w2", "Factor3_w4",
                                        "Factor4_w0", "Factor4_w1", "Factor4_w2", "Factor4_w4",
                                        "Factor5_w0", "Factor5_w1", "Factor5_w2", "Factor5_w4"),
                         variable.name ="TimeAfterInjury", 
                         value.name ="FactorsValues")

Factors_long <- separate(data = Factors_long, col = TimeAfterInjury, 
                             into = c("Factors", "TimeAfterInjury"), sep = "\\_")

library(party)
UCP_Factors <- ctree(TLAMS52 ~ AGE + SEXCD + ASIMPC01_A + SPLVL1 + TimeAfterInjury + FactorsValues, data = na.omit(Factors_long))
plot(UCP_Factors)

#UCP is not working for long format so I used the wide format for each timepoint
#Week 0
Factor1 <- cfa_week0[, c("RBC00", "HGB00", "HCT00")]
Factors_W0 <- do.call("cbind", list(MS52_raw, p_week0, Other_Variables, Factor1))

#Coding for death outcome
Factors_W0$DEATHRPD[is.na(Factors_W0$DEATHRPD)] <- 0
Factors_W0$DEATHRPD[Factors_W0$DEATHRPD == 2] <- 1

#Having to code for death first before transforing it into factors values

#Turning variables into factor
str(Factors_W0)
Factors_W0 <- Factors_W0 %>% mutate_if(is.character,as.factor) #They only take factor values
Factors_W0$SEXCD <- as.factor(Factors_W0$SEXCD) 
Factors_W0$DEATHRPD <- as.factor(Factors_W0$DEATHRPD)

#ctree can't take missing values in response!
Factors_W0 <- subset(Factors_W0, !is.na(TLAMS52))

#UCP for LMS52
UCP_Factors_W0 <- ctree(TLAMS52 ~ Factor1_w0 + Factor2_w0 + SEXCD + AGE + SPLVL1 + ASIMPC01_A +MAP + HEART + RESPRT + TEMPC + BOLDOSE1 + INFDOSE1, data = Factors_W0)
plot(UCP_Factors_W0)

#GLM for death outcome
x <- glm(DEATHRPD ~ Factor1_w0 + Factor2_w0 + SEXCD + AGE + SPLVL1 + ASIMPC01_A + MAP + HEART + RESPRT + TEMPC, data = Factors_W0, family = "binomial")
summary(x)
anova(x, test = "Chisq")

#Subset ASI-A grade to see death outcome
Factors_W0_A <- Factors_W0[Factors_W0$ASIMPC01_A %in% "A", ]

x <- glm(DEATHRPD ~ Factor1_w0 + Factor2_w0 + SEXCD + AGE + SPLVL1, data = Factors_W0_A, family = "binomial")
summary(x)
anova(x, test = "Chisq")

UCP_Factors_W0_A <- ctree(DEATHRPD ~ Factor1_w0 + Factor2_w0 + SEXCD + AGE + SPLVL1, data = Factors_W0_A)
plot(UCP_Factors_W0_A)

#UCP for week 4
Factors_W4 <- do.call("cbind", list(MS52_raw, p_week4, Other_Variables, Factor1))

str(Factors_W4)
Factors_W4 <- Factors_W4 %>% mutate_if(is.character,as.factor) #They only take factor values
Factors_W4$SEXCD <- as.factor(Factors_W4$SEXCD) 

Factors_W4 <- subset(Factors_W4, !is.na(TLAMS52))

UCP_Factors_W4 <- ctree(TLAMS52 ~ Factor1 + Factor2 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_W4)
plot(UCP_Factors_W4)


#UCP for week1-4
Factors_Overtime <- Factors_Overtime %>% mutate_if(is.character,as.factor)
Factors_Overtime$SEXCD <- as.factor(Factors_Overtime$SEXCD) 

str(Factors_Overtime)

Factors_Overtime <- subset(Factors_Overtime, !is.na(TLAMS52))

UCP_Factors_Overtime <- partykit::ctree(TLAMS52 ~ Factor1_w0 + Factor1_w1 + Factor1_w2 + Factor1_w4 + SEXCD + SPLVL1 + AGE + ASIMPC01_A, data = Factors_Overtime)
plot(UCP_Factors_Overtime)


library(strucchange)
sctest(UCP_Factors_Overtime, node = 1)

#Growth Curve Latent 

Model5 <- '
    i =~ 1*Factor1_w0 + 1*Factor1_w1 + 1*Factor1_w2 + 1*Factor1_w4
    s =~ 0*Factor1_w0 + 1*Factor1_w1 + 2*Factor1_w2 + 3*Factor1_w4
    '

curve <- growth(Model5, Factors_Overtime)
summary(curve)

my_curve <- lavaan(Model5, Factors_Overtime)
summary(my_curve)

tree <- semtree(curve, Factors_Overtime, control=semtree.control(alpha = 0.05, bonferroni = TRUE, exclude.heywood = TRUE))
plot(tree)

#scaling data of week 8
cfa1_week8 <- scale(cfa_week8, center = TRUE, scale = TRUE)
cfa1_week8 <- data.frame(cfa1_week8)
cfa1_week8

scree <- fa.parallel(cfa1_week8, fm = "ml", fa = "fa")
scree

fa <- fa(cfa1_week8, nfactors=6, rotate="varimax", fm="ml")
fa.diagram(fa)
