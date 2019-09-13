## Change value of base to path to data and scripts
base <- "/Users/benbarnes2/Documents/DGEpi/Modelling Workshop/Workshop APC/R code and data"

if(!exists("finBladPop")) {
  source(file.path(base, "APC data prep.R"))
}
p0 <- 1976
c0 <- 1906

## Plot rates by calendar year and period, group by age group
## Bladder cancer
ggplot(finBladMenAPC[meanAge >= 30 & perBeg > 1955], aes(x = cohMid, y = newCrudeRate,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

ggplot(finBladMenAPC[meanAge >= 30 & perBeg > 1955], aes(x = perBeg, y = newCrudeRate,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

## Lung cancer among women
ggplot(finLungWomenAPC[meanAge >= 50 & perBeg > 1955], aes(x = cohMid, y = newCrudeRate,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

ggplot(finLungWomenAPC[meanAge >= 50 & perBeg > 1955], aes(x = perBeg, y = newCrudeRate,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

## Lung cancer among men
ggplot(finLungMenAPC[meanAge >= 50 & perBeg > 1955], aes(x = cohMid, y = newCrudeRate,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

ggplot(finLungMenAPC[meanAge >= 50 & perBeg > 1955], aes(x = perBeg, y = newCrudeRate,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")


## Create some models
## Age-cohort factor model with 5-year period data
glmBlad <- glm(cases ~ ageGroup + cohMid1925 + offset(log(pop)) - 1,
  data = finBladMenAPC[meanAge >= 30 & meanAge < 80], family = poisson(link = "log"))
summary(glmBlad)
exp(cbind(coef(glmBlad), confint(glmBlad)))
## Predictions
finBladMenAPC[meanAge >= 30 & meanAge < 80,
  predInc := predict(glmBlad, newdata = finBladMenAPC[meanAge >= 30 & meanAge < 80],
    type = "response") / pop * 1e5]

ggplot(finBladMenAPC[meanAge >= 30 & meanAge < 80], aes(x = cohMid, y = predInc,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

## Age-period model using splines and 1-year periods
## glm function
agePerGLM <- glm(D ~ splines::bs(A, df = 8) + splines::bs(P - p0, df = 8) +
  offset(log(Y)),
  data = finBladPopMen, family = poisson(link = "log"))
summary(agePerGLM)

finBladPopMen[, predAPglm := predict(agePerGLM,
  newdata = finBladPopMen, type = "response")]

ggplot(finBladPopMen[meanAge >= 30 & meanAge < 80], aes(x = P, y = predAPglm,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

ggplot(finBladPopMen[meanAge >= 30 & meanAge < 80], aes(x = cohMid, y = predAPglm,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

## apc.fit function
agePeriod <- apc.fit(finBladPopMen[meanAge >= 30 & meanAge < 80], dist = "poisson",
  model = "bs", npar = c(A = 9, P = 9, C = 8), parm = "AP-C",
  scale = 1e5, ref.p = p0)

## This is the age-period model component:
summary(agePeriod[["Model"]][[1]])

## Cohort effect estimated based on residuals from age-period model
plot(agePeriod, lty = c(1, 1, 3), col = "red", shade = TRUE)

finBladPopMen[meanAge >= 30 & meanAge < 80,
  predAPepi := predict(agePeriod[["Model"]][[1]], type = "response")]
ggplot(finBladPopMen[meanAge >= 30 & meanAge < 80], aes(x = P, y = predAPepi,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

ggplot(finBladPopMen[meanAge >= 30 & meanAge < 80], aes(x = cohMid, y = predAPepi,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

## Age-Cohort model
ageCohort <- apc.fit(finBladPopMen[meanAge >= 30 & meanAge < 80], dist = "poisson",
  model = "bs", npar = c(A = 9, P = 9, C = 8), parm = "AC-P",
  scale = 1e5)

plot(ageCohort, lty = c(1, 1, 3), shade = TRUE)

## Compare drift parameters (they are the same)
agePeriod$Drift
ageCohort$Drift

## age-period-cohort model

agePerCoh <- apc.fit(finBladPopMen[meanAge >= 30 & meanAge < 80], dist = "poisson",
  model = "bs", npar = c(A = 9, P = 9, C = 8), parm = "APC",
  scale = 1e5, ref.p = 1982.5)
plot(agePerCoh, shade = TRUE)
agePerCoh[["Age"]]
agePerCoh[["Per"]]
agePerCoh[["Coh"]]

## age-cohort-period model

ageCohPer <- apc.fit(finBladPopMen[meanAge >= 30 & meanAge < 80], dist = "poisson",
  model = "bs", npar = c(A = 9, P = 9, C = 8), parm = "ACP",
  scale = 1e5)
plot(ageCohPer, shade = TRUE)

## compare fitted values from APC and ACP models
finBladPopMen[meanAge >= 30 & meanAge < 80, fittedAPC := predict(agePerCoh[["Model"]],
  type = "response") / pop * 1e5]
finBladPopMen[meanAge >= 30 & meanAge < 80, fittedACP := predict(ageCohPer[["Model"]],
  type = "response") / pop * 1e5]
## fitted values are the same
finBladPopMen[meanAge >= 30 & meanAge < 80,
  list(ageGroup, period5, cohort5, fittedAPC, fittedACP, newCrudeRate)]

identical(agePerCoh[["Model"]]$model$D,
  finBladPopMen[meanAge >= 30 & meanAge < 80, D])

