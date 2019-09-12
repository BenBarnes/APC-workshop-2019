## Projections

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

