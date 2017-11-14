efa <-  function(x) {
  scree <- fa.parallel(x[-1,-2, -3], fm = "ml", fa = "fa")
  fa <- fa.diagram(fa(x[-1, -2, -3], nfactors=6, rotate="oblimin", fm="ml"))
  return(list(scree, fa))
}

efa(cfa1_week0)

library(psych)
scree <- fa.parallel(cfa1_week0[-1:-3], fm = "ml", fa = "fa")
scree


f1 <- lm(Outcome ~ Factor1, data = p_week0)
f2 <- lm(Outcome ~ Factor2, data = p_week0)
f3 <- lm(Outcome ~ Factor3, data = p_week0)
f4 <- lm(Outcome ~ Factor4, data = p_week0)
f5 <- lm(Outcome ~ Factor5, data = p_week0)

factor_lm <- function(x) {
  f1 <- summary(lm(Outcome ~ Factor1, data = x))
  f2 <- summary(lm(Outcome ~ Factor2, data = x))
  f3 <- summary(lm(Outcome ~ Factor3, data = x))
  f4 <- summary(lm(Outcome ~ Factor4, data = x))
  f5 <- summary(lm(Outcome ~ Factor5, data = x))
  return(list(f1,f2,f3,f4,f5))
}

factor_lm(p_week0)

cors <- cor(subset(mtcars, select = c(-gear,-carb)))
cors2 <- sapply(subset(mtcars, select = c(-gear,-carb)), cor(mtcars,use = "pairwise.complete.obs"), y=mtcars$mpg)
cors2


sapply(cfa2_week0, cor, y=cfa2_week0$Outcome)

cfa2_week0 <- merge(cfa1_week0, p_week0, how = 'outer')
cfa2_week0

one <- subset(cfa2_week0, select = c(clusters))

w <- cbind(one, p_week0)

ab <- rbind.fill(cfa2_week0, p_week0)

cfa2_week0 <- cbind(cfa1_week0, p_week0)
cfa2_week0

c <- cor(factors, use="pairwise.complete.obs")
c <- as.data.frame(c)
b <- c[1:3]
