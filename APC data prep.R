library(data.table)
library(Epi)
library(ggplot2)
library(XML)

## Change value of base to path to data and scripts
if(!exists("base", envir = .GlobalEnv)) {
  base <- "/Users/benbarnes2/Documents/DGEpi/Modelling Workshop/Workshop APC/R code and data"
}

theAges <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
  "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85plus")
theSexes <- c("men", "women")
theDiagnoses <- c(Lung = 180, Bladder = 300, Prostate = 261)

## Finnish lung cancer data among women:
getFinData <- function(diagnosis, sex) {
  theDiag <- theDiagnoses[match(diagnosis, names(theDiagnoses))]
  res <- lapply(seq_len(18), function(theAge) {# browser()
    theURL <- sprintf(
      "http://www-dep.iarc.fr/NORDCAN/english/Table1.asp?cancer=%i&registry=246&sYear=1953&eYear=2016&sex=%i&type=0&age_from=%i&age_to=%i&submit=Execute",
      theDiag, sex, theAge, theAge)
    theTabs <- getNodeSet(htmlParse(theURL), "//table")
    resTab <- readHTMLTable(theTabs[[2]], stringsAsFactors = FALSE)
    setDT(resTab)
    resTab[, ageGroup := theAges[theAge]]
  })
  res2 <- rbindlist(res)[!grepl("-", Year)]
  res2[, Sex := theSexes[sex]]
  res2[, Diagnosis := diagnosis]
  res2[, Numbers := as.numeric(Numbers)]
  res2[, Year := as.numeric(Year)]
  return(res2)
}

## Read in Finnish cancer data
lungIncWomen <- getFinData(diagnosis = "Lung", sex = 2)
bladIncMen <- getFinData(diagnosis = "Bladder", sex = 1)
lungIncMen <- getFinData(diagnosis = "Lung", sex = 1)

## Read in Finnish male population data
finPop <- fread(file.path(base, "Finland population.csv"), skip = 2)
## Melt to long format
mFinPop <- melt(finPop, id.vars = "Year", value.name = "pop")
## Determine sex and age group values
mFinPop[, Sex := ifelse(grepl("Males", variable), "men", "women")]
mFinPop[, ageGroup := sub(".+ ([0-9]{1,2}) - ([0-9]{1,2}).+", "\\1-\\2", variable)]
mFinPop[grepl("85", ageGroup), ageGroup := "85plus"]
mFinPop[grepl("- 4", ageGroup), ageGroup := "0-4"]

## Calculate mid year of 5-year cohort based on calendar year and age group
mFinPop[, meanAge := as.numeric(sub("^([0-9]{1,2})[-p].+", "\\1", ageGroup)) + 3]
mFinPop[, cohMid := Year - meanAge]
mFinPop[, cohMid1925 := relevel(factor(cohMid), ref = "1925")]

## Create 5-year cohorts and 5-year periods add up incidence data
mFinPop[, period5 := cut(Year, seq(min(Year), max(Year), by = 5), right = FALSE)]
mFinPop[, cohort5 := cut(cohMid, seq(min(cohMid), max(cohMid), by = 5),
  right = FALSE)]
mFinPop[, cohort5_1925 := relevel(cohort5, ref = "[1925,1930)")]
stopifnot(all(mFinPop[, length(unique(ageGroup)) == 1,
  by = list(period5, cohort5)]$V1))

## Merge incidence and population data
finBladPopMen <- merge(mFinPop, bladIncMen, by.x = c("Year", "Sex", "ageGroup"),
  by.y = c("Year", "Sex", "ageGroup"))
## Re-calculate age-specific rates
finBladPopMen[, newCrudeRate := Numbers / pop * 1e5]

finLungPopWomen <- merge(mFinPop, lungIncWomen, by.x = c("Year", "Sex", "ageGroup"),
  by.y = c("Year", "Sex", "ageGroup"))

finBladMenAPC <- finBladPopMen[!is.na(period5), list(cases = sum(Numbers), pop = sum(pop),
  perBeg = min(Year), cohMid = max(cohMid)),
  keyby = list(ageGroup, meanAge, period5, cohort5)]
finLungWomenAPC <- finLungPopWomen[!is.na(period5), list(cases = sum(Numbers), pop = sum(pop),
  perBeg = min(Year), cohMid = max(cohMid)),
  keyby = list(ageGroup, meanAge, period5, cohort5)]
## Re-calculate age-specific rates
finBladMenAPC[, newCrudeRate := cases / pop * 1e5]
finLungWomenAPC[, newCrudeRate := cases / pop * 1e5]
## set reference level of cohorts to cohort ending in 1928
finBladMenAPC[, cohMid1925 := relevel(factor(cohort5), ref = "[1925,1930)")]
finLungWomenAPC[, cohMid1925 := relevel(factor(cohort5), ref = "[1925,1930)")]

## Variables required for apc.fit function
finBladPopMen[, c("A", "P", "D", "Y") := list(meanAge, Year, Numbers, pop)]
finLungPopWomen[, c("A", "P", "D", "Y") := list(meanAge, Year, Numbers, pop)]
finBladMenAPC[, c("A", "P", "D", "Y") := list(meanAge, perBeg + 3, cases, pop)]
finLungWomenAPC[, c("A", "P", "D", "Y") := list(meanAge, perBeg + 3, cases, pop)]

