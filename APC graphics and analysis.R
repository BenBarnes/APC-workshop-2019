## Change value of base to path to data and scripts
base <- "/Users/benbarnes2/Documents/DGEpi/Modelling Workshop/Workshop APC/R code and data"

if(!exists("finBladPop")) {
  source(file.path(base, "APC data prep.R"))
}

## Plot rates by calendar year and period, group by age group
## Bladder cancer
ggplot(finBladMenAPC[meanAge >= 50 & perBeg > 1955], aes(x = cohMid, y = newCrudeRate,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

ggplot(finBladMenAPC[meanAge >= 50 & perBeg > 1955], aes(x = perBeg, y = newCrudeRate,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

## Lung cancer among women
ggplot(finLungWomenAPC[meanAge >= 50 & perBeg > 1955], aes(x = cohMid, y = newCrudeRate,
  color = ageGroup)) +
  geom_line() + geom_point() + coord_trans(y = "log2")

ggplot(finLungWomenAPC[meanAge >= 50 & perBeg > 1955], aes(x = perBeg, y = newCrudeRate,
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

agePeriod <- apc.fit(finBladPopMen[meanAge >= 30 & meanAge < 80], dist = "poisson",
  model = "bs", npar = c(A = 9, P = 9, C = 8), parm = "AP-C",
  scale = 1e5, ref.p = 1982.5)

## This is the age-period model component:
summary(agePeriod[["Model"]][[1]])

## Cohort effect estimated based on residuals from age-period model
plot(agePeriod, lty = c(1, 1, 3), col = "red", shade = TRUE)

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

## WORK IN PROGRESS
## Projections with new data
## This will only be applicable to age-period models
p0 <- 1976
c0 <- 1906
finBladTrain <- finBladPopMen[Year <= 1985]

ap_c_vars <- function(theData, ref.p, ref.c, returnWhat = "all",
  knotsA, knotsP, knotsC, bKnotsA, bKnotsP, bKnotsC) {
  med <- function(x, y) {
    o <- order(x)
    a <- y[o]
    names(a) <- x[o]
    return(as.numeric(names(a[cumsum(a)/sum(a) > 0.5][1])))
  }
  has.pref <- !missing(ref.p)
  D <- theData$D
  P <- theData$P
  A <- theData$A
  p0 <- ifelse(has.pref, ref.p, med(P, D))
  has.cref <- !missing(ref.c)
  c0 <- ifelse(has.cref, ref.c, med(P - A, D))
  Y <- theData$Y
  if(!missing(knotsA)) {
    MA0 <- theData[, splines::bs(A, knots = knotsA, Boundary.knots = bKnotsA, degree = 3)]
  } else {
    MA0 <- theData[, splines::bs(A, df = 9, degree = 3)]
  }
  MA <- cbind(1, MA0)
  if(!missing(knotsP)) {
    MP <- theData[, splines::bs(P, knots = knotsP, Boundary.knots = bKnotsP, degree = 3)]
  } else {
    MP <- theData[, splines::bs(P, df = 9, degree = 3)]
  }
  if(!missing(knotsC)) {
    MC <- theData[, splines::bs(P - A, knots = knotsC, Boundary.knots = bKnotsC, degree = 3)]
  } else {
    MC <- theData[, splines::bs(P - A, df = 8, degree = 3)]
  }
  Rp <- splines::bs(p0, knots = attr(MP, "knots"), Boundary.knots = attr(MP, 
    "Boundary.knots"), degree = attr(MP, "degree"))
  Rc <- splines::bs(c0, knots = attr(MC, "knots"), Boundary.knots = attr(MC, 
    "Boundary.knots"), degree = attr(MC, "degree"))
  wt <- Y
  xP <- Epi::detrend(rbind(Rp, MP), c(p0, P), weight = c(0, wt))
  MPr <- xP[-1, , drop = FALSE] - has.pref * xP[rep(1, nrow(MP)), 
    , drop = FALSE]
  lP <- cbind(P - p0, MPr)
  if(returnWhat == "all"){
    out <- list(D = D, MA = MA, P = P, p0 = p0, MP = MP, MC = MC, Y = Y, lP = lP,
      knotsA = attr(MA0, "knots"), knotsP = attr(MP, "knots"), knotsC = attr(MC, "knots"),
      bKnotsA = attr(MA0, "Boundary.knots"), bKnotsP = attr(MP, "Boundary.knots"),
      bKnotsC = attr(MC, "Boundary.knots"))
  } else {
    if(returnWhat == "MA") {
      out <- MA
    } else {
      if(returnWhat == "lP") {
        out <- lP
      }
    }
  }
  return(out)
}

ap_c_mods <- function(trainVars) {
  D <- trainVars$D
  MA <- trainVars$MA
  P <- trainVars$P
  p0 <- trainVars$p0
  MP <- trainVars$MP
  MC <- trainVars$MC
  Y <- trainVars$Y
  lP <- trainVars$lP
  m.APC <- glm(D ~ MA + I(P - p0) + MP + MC + offset(log(Y)), family = poisson)
  m.AP <- update(m.APC, . ~ . - MC)
  m.0 <- update(m.AP, . ~ . - MA - I(P - p0) - MP)
  ap <- update(m.0, . ~ . - 1 + MA + lP)
  return(ap)
}

agePerTrain <- apc.fit(finBladTrain, dist = "poisson", model = "bs",
  npar = c(A = 9, P = 9, C = 8), parm = "AP-C",
  scale = 1e5, ref.p = p0)

trainVars <- ap_c_vars(theData = finBladTrain, ref.p = p0)
trainRepro <- ap_c_mods(trainVars)
identical(coef(agePerTrain[["Model"]][[1]]), coef(trainRepro))

## Have to deal with offset
# Remove offset from formula and terms
modFormEnv <- attr(trainRepro$formula, ".Environment")
trainRepro$formula <- as.formula("D ~ MA + lP - 1",
  env = modFormEnv)
trainRepro$terms <- terms(trainRepro$formula)


## Check to see whether link predictions are the same with altered environments
myPred <- predict(trainRepro, type = "link")
origPred <- predict(agePerTrain$Model[[1]], type = "link")
head(cbind(myPred, origPred))

## Try to assign spline bases from full dataset to model environment
fullVars <- ap_c_vars(finBladPopMen, ref.p = p0, ref.c = c0, returnWhat = "all",
  knotsA = trainVars$knotsA, knotsP = trainVars$knotsP, knotsC = trainVars$knotsC,
  bKnotsA = trainVars$bKnotsA, bKnotsP = trainVars$knotsP, bKnotsC = trainVars$knotsC)
invisible(lapply(seq_along(fullVars), function(theVar) {# browser()
  assign(names(fullVars)[theVar], fullVars[[theVar]], modFormEnv)
}))

length(get("lP", modFormEnv))

length(bladPred <- predict(trainRepro, newdata = finBladPop, type = "response"))

finBladPop[, pred1985 := predict(trainRepro, newdata = finBladPop,
  type = "response") * pop]
finBladPop[Year <= 1985, predTrain := predict(agePerTrain[["Model"]][[1]],
  type = "response")]

finBladPop[Year <= 1985]

# data <- .External2(C_modelframe, myForm, rownames, variables[1:2], 
#   varnames[1:2], extras, extranames, subset, na.action)

