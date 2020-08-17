# Maximum-Likelihood estimation of R
### Michelle Kendall and Christophe Fraser



This code performs a Maximum-Likelihood estimation of R over fixed time periods before and after the Test and Trace programme, as reported in our [preprint](https://www.medrxiv.org/content/10.1101/2020.07.12.20151753v1.article-info).

## Setup

### Load packages


```r
library(EpiEstim)
library(incidence)
library(ggplot2)
library(stats4)
library(tidyverse)
library(RColorBrewer)
```

## Setup

### Load data

From the working directory, fetch the daily case data:


```r
dat.UK.p1 <- read_csv("data/pillar1_case_data.csv") # local area pillar 1 data
dat.UK.p2 <- read.csv("data/pillar2_case_data.csv") # national pillars data including pillar 2
dat.IoW.p2 <- read.csv("data/pillar2_case_data_IoW.csv") # Isle of Wight pillar 2 data
dat.UK.both <- read.csv("data/combined_pillars_case_data.csv") # local combined pillars data
```

and the population size data:


```r
population.data <- read.csv("data/population_by_region.csv", stringsAsFactors = FALSE)
```

### Cleaning data

Extract data for England and Isle of Wight (Pillar 1)


```r
# Extract relevant data:
dat.Eng.p1 <- dat.UK.p1 %>% filter(Country == "England") # England only
dat.IoW.p1 <- dat.Eng.p1 %>% filter(Area == "Isle of Wight") # IoW only

dat.Eng.both <- dat.UK.both %>% filter(`Area.type` == "Upper tier local authority") # UTLAs only
dat.Eng.both.clean <- dat.Eng.both %>% select("Date"=`Specimen.date`, "Area"=`Area.name`, "Cases"=`Daily.lab.confirmed.cases`) # only the necessary columns
dat.Eng.both.clean$Date <- as.Date(dat.Eng.both.clean$Date) # restore "date" property
```


## Analysis

Calculate incidence for each area and overall:


```r
areas.alphabetical <- sort(unique(dat.Eng.both.clean$Area)) # get a list of utlas in alphabetical order
IoW.index <- which(areas.alphabetical == "Isle of Wight")

earliest.date <- min(as.Date(dat.Eng.both.clean$Date))
last.date <- max(as.Date(dat.Eng.both.clean$Date))
all.dates <- seq(earliest.date, last.date, by=1)

# Get incidence by area, filling in any missing dates with 0 cases
incidence.each.area <- sapply(areas.alphabetical, function(area) {
  dat.area <- dat.Eng.both.clean %>% filter(Area == area)
  
  dat.area <- rbind.data.frame(dat.area,
    cbind.data.frame("Date" = as.Date(setdiff(all.dates, dat.area$Date), origin = "1970-01-01"),
                     "Area" = area,
                     "Cases" = 0)
  )

  dat.area$Cases[order(dat.area$Date)]
  
})

# extract the incidence table for the Isle of Wight
IoW.incidence <- cbind.data.frame(
  "Date" = all.dates,
  "Cases" = incidence.each.area[,IoW.index]
)

# and overall for England:
total.incidence <- cbind.data.frame(
  "Date" = all.dates,
  "Cases" = sapply(1:nrow(incidence.each.area), function(d) sum(incidence.each.area[d,]))
)

# these can be visualised for example like this:
TotalIncidencePlot <- ggplot(total.incidence, aes(x=Date, y=Cases)) +
  geom_col(fill= brewer.pal(6,"Set2")[[3]]) +
  scale_x_date(date_breaks = "2 weeks" , date_labels = "%b %d", limits=c(as.Date("2020-04-01"), last.date - 7)) +
  ylab("Daily confirmed COVID-19\ncases by swab date") +
  theme_bw(base_size = 16) 

TotalIncidencePlot
```

![](figure/analysis.both.1-1.png)

For England and for the Isle of Wight create their incidence line lists and `incidence` objects using data from 5 days after the TT launch and ending 5 days before the end (to allow for missing data as cases where there are delays in reporting cases.)


```r
start.date.IoW.both <- as.Date("2020-05-10")
start.date.Eng.both <- as.Date("2020-05-23")
end.date.both <- last.date - 5 

IoW.both.linelist <- IoW.incidence$Date[1]
for(row in 1:nrow(IoW.incidence)){
  if(IoW.incidence$Cases[row]>0){
    for(case in 1:IoW.incidence$Cases[row]){
      IoW.both.linelist <- c(IoW.both.linelist, IoW.incidence$Date[row])
    }
  }
}
IoW.both.linelist <- IoW.both.linelist[-1]

total.incidence$Cases[which(is.na(total.incidence$Cases))] <- 0
Eng.both.linelist <- total.incidence$Date[1]
for(row in 1:nrow(total.incidence)){
  if(total.incidence$Cases[row]>0){
    Eng.both.linelist <- c(Eng.both.linelist, rep(total.incidence$Date[row], total.incidence$Cases[row]))
  }
}
Eng.both.linelist <- Eng.both.linelist[-1]

IoW.both.incidence <- incidence(IoW.both.linelist, 
                              first_date = start.date.IoW.both,
                              last_date = end.date.both,
                              standard = FALSE)

Eng.both.incidence <- incidence(Eng.both.linelist, 
                             first_date = start.date.Eng.both,
                             last_date = end.date.both,
                             standard = FALSE)
```

Prepare log-likelihood function:

```r
minusloglik <- function(i0, r, start.date, incidence.dates, incidence.counts) {
  ll <- 0
  for(row in 1:length(incidence.dates)){
    model <- i0*exp(-r*as.numeric(difftime(start.date,incidence.dates[row],units=c("days"))))
    ll <- ll + model - incidence.counts[row]*log(model)
  }
  return(ll)
}
```

Get Maximum-Likelihood estimates:


```r
fit.mle.IoW.both <- mle(minusloglik, 
                      start=list(i0=10, r=-0.01), 
                      fixed=list(start.date=start.date.IoW.both,
                                 incidence.dates=IoW.both.incidence$dates,
                                 incidence.counts=IoW.both.incidence$counts)
)
summary(fit.mle.IoW.both)
```

```
## Maximum likelihood estimation
## 
## Call:
## mle(minuslogl = minusloglik, start = list(i0 = 10, r = -0.01), 
##     fixed = list(start.date = start.date.IoW.both, incidence.dates = IoW.both.incidence$dates, 
##         incidence.counts = IoW.both.incidence$counts))
## 
## Coefficients:
##       Estimate  Std. Error
## i0  9.62302144 1.236653064
## r  -0.07209213 0.008048002
## 
## -2 log L: -110.3964
```

```r
fit.mle.Eng.both <- mle(minusloglik,
                     start=list(i0=1420, r=-0.02), 
                     fixed=list(start.date=start.date.Eng.both,
                                incidence.dates=Eng.both.incidence$dates,
                                incidence.counts=Eng.both.incidence$counts))
summary(fit.mle.Eng.both)
```

```
## Maximum likelihood estimation
## 
## Call:
## mle(minuslogl = minusloglik, start = list(i0 = 1420, r = -0.02), 
##     fixed = list(start.date = start.date.Eng.both, incidence.dates = Eng.both.incidence$dates, 
##         incidence.counts = Eng.both.incidence$counts))
## 
## Coefficients:
##         Estimate  Std. Error
## i0 1420.15949911 1.35301e+01
## r    -0.02283999 5.38094e-04
## 
## -2 log L: -410782.2
```


Similarly for Pillar 1 alone:

```r
dat.IoW.p1$Incidence <- rep(0,nrow(dat.IoW.p1)) 
dat.IoW.p1$Incidence[1] <- dat.IoW.p1$TotalCases[1]
for(row in 2:nrow(dat.IoW.p1)) {
  dat.IoW.p1$Incidence[row] <- dat.IoW.p1$TotalCases[row] - dat.IoW.p1$TotalCases[row-1]
}


dat.Eng.p1 <- dat.UK.p1 %>% filter(Country == "England")

incidence.all.areas.p1 <- sapply(areas.alphabetical, function(area) {
  dat.area <- dat.Eng.p1 %>% filter(Area == area)
  
  new.dat.area <- data.frame("Date" = unique(dat.Eng.p1$Date))
  new.dat.area$TotalCases <- sapply(new.dat.area$Date, function(d) (dat.area %>% filter(Date == d))$TotalCases)
  # fill in the blanks
  if (length(new.dat.area$TotalCases[[1]])==0) new.dat.area$TotalCases[[1]] <- 0
  for (i in 1:nrow(new.dat.area)) {
    if (length(new.dat.area$TotalCases[[i]])==0) new.dat.area$TotalCases[[i]] <- new.dat.area$TotalCases[[i-1]]
  }
  
  new.dat.area$Incidence <- rep(0,nrow(new.dat.area))
  new.dat.area$Incidence[[1]] <- new.dat.area$TotalCases[[1]]
  for(row in 2:nrow(new.dat.area)) new.dat.area$Incidence[[row]] <- new.dat.area$TotalCases[[row]] - new.dat.area$TotalCases[[row-1]]
  
  new.dat.area$Incidence
})

Eng.Pillar1 <- cbind.data.frame(
  "Date" = unique(dat.Eng.p1$Date),
  "Incidence" = sapply(1:nrow(incidence.all.areas.p1), function(row) {
    sum(incidence.all.areas.p1[row,], na.rm=TRUE)
  })
)


# linelists
IoW.p1.linelist <- dat.IoW.p1$Date[1]
for(row in 1:nrow(dat.IoW.p1)){
  if(dat.IoW.p1$Incidence[row]>0){
    for(case in 1:dat.IoW.p1$Incidence[row]){
      IoW.p1.linelist <- c(IoW.p1.linelist, dat.IoW.p1$Date[row])
    }
  }
}
IoW.p1.linelist <- IoW.p1.linelist[-1]

Eng.p1.linelist <- Eng.Pillar1$Date[1]
for(row in 1:nrow(Eng.Pillar1)){
  if ((!is.na(Eng.Pillar1$Incidence[row]))&&(Eng.Pillar1$Incidence[row]>0)){
    Eng.p1.linelist <- c(Eng.p1.linelist, rep(Eng.Pillar1$Date[row], Eng.Pillar1$Incidence[row]))
  }
}
Eng.p1.linelist <- Eng.p1.linelist[-1]
```

And Pillar 2:

```r
dat.IoW.p2$specimen_date <- as.Date(as.character(dat.IoW.p2$specimen_date), "%d-%b-%y")
dat.IoW.p2$Pillar2 <- sapply(dat.IoW.p2$Pillar2, function(x)(if(!is.na(x)) return(x) else return(0)))

dat.UK.p2$Date.of.activity <- as.Date(dat.UK.p2$Date.of.activity, "%d/%m/%y")
UKPillar2 <- dat.UK.p2 %>% filter(Pillar == "Pillar 2 (excluding Wales)")
stopifnot(unique(UKPillar2$Nation) == "UK") # check we've just got the UK entries
UKPillar2$Daily.number.of.positive.cases <- sapply(UKPillar2$Daily.number.of.positive.cases, function(x)(if(!is.na(x)) return(x) else return(0)))

IoW.p2.linelist <- dat.IoW.p2$specimen_date[1]
for(row in 1:nrow(dat.IoW.p2)){
  if(dat.IoW.p2$Pillar2[row]>0){
    for(case in 1:dat.IoW.p2$Pillar2[row]){
      IoW.p2.linelist <- c(IoW.p2.linelist, dat.IoW.p2$specimen_date[row])
    }
  }
}
IoW.p2.linelist <- IoW.p2.linelist[-1]

UK.p2.linelist <- UKPillar2$Date.of.activity[1]
for(row in 1:nrow(UKPillar2)){
  if(UKPillar2$Daily.number.of.positive.cases[row]>0){
    # this is faster:
    UK.p2.linelist <- c(UK.p2.linelist, rep(UKPillar2$Date.of.activity[row], UKPillar2$Daily.number.of.positive.cases[row]))
    # for(case in 1:UKPillar2$Daily.number.of.positive.cases[row]){
    #   UK.linelist <- c(UK.linelist, UKPillar2$Date.of.activity[row])
    # }
  }
}
UK.p2.linelist <- UK.p2.linelist[-1]
```


Get `incidence` objects and Maximum-Likelihood estimates for the individual pillars data, starting 10 days after TT launch for Pillar 1 and 5 days after for Pillar 2 (to reflect the different typical times between symptom onset and swab).


```r
start.date.IoW.p1 <- as.Date("2020-05-15")
start.date.IoW.p2 <- as.Date("2020-05-10")
start.date.Eng.p1 <- as.Date("2020-05-28")
start.date.UK.p2 <- as.Date("2020-05-23")

end.date.IoW.p1 <- max(Eng.Pillar1$Date) - 5 
end.date.Eng.p1 <- max(Eng.Pillar1$Date) - 5
end.date.UK.p2 <- max(UKPillar2$Date.of.activity, na.rm = TRUE) - 5 
end.date.IoW.p2 <- end.date.UK.p2 # there were no more cases after the 13th

IoW.p1.incidence <- incidence(IoW.p1.linelist, 
                              first_date = start.date.IoW.p1,
                              last_date = end.date.IoW.p1,
                              standard = FALSE)

Eng.p1.incidence <- incidence(Eng.p1.linelist, 
                             first_date = start.date.Eng.p1,
                             last_date = end.date.Eng.p1,
                             standard = FALSE)

IoW.p2.incidence <- incidence(IoW.p2.linelist, 
                           first_date = start.date.IoW.p2,
                           last_date = end.date.IoW.p2,
                           standard = FALSE)

UK.p2.incidence <- incidence(UK.p2.linelist, 
                          first_date = start.date.UK.p2,
                          last_date = end.date.UK.p2,
                          standard = FALSE)
```


```r
fit.mle.IoW.p1 <- mle(minusloglik, 
                      start=list(i0=10, r=-0.01), 
                      fixed=list(start.date=start.date.IoW.p1,
                                 incidence.dates=IoW.p1.incidence$dates,
                                 incidence.counts=IoW.p1.incidence$counts)
)
summary(fit.mle.IoW.p1)
```

```
## Maximum likelihood estimation
## 
## Call:
## mle(minuslogl = minusloglik, start = list(i0 = 10, r = -0.01), 
##     fixed = list(start.date = start.date.IoW.p1, incidence.dates = IoW.p1.incidence$dates, 
##         incidence.counts = IoW.p1.incidence$counts))
## 
## Coefficients:
##      Estimate Std. Error
## i0  3.6521234 0.92358594
## r  -0.1245396 0.02442024
## 
## -2 log L: 37.98268
```

```r
fit.mle.Eng.p1 <- mle(minusloglik,
                     start=list(i0=10, r=-0.03), 
                     fixed=list(start.date=start.date.Eng.p1,
                                incidence.dates=Eng.p1.incidence$dates,
                                incidence.counts=Eng.p1.incidence$counts))
summary(fit.mle.Eng.p1)
```

```
## Maximum likelihood estimation
## 
## Call:
## mle(minuslogl = minusloglik, start = list(i0 = 10, r = -0.03), 
##     fixed = list(start.date = start.date.Eng.p1, incidence.dates = Eng.p1.incidence$dates, 
##         incidence.counts = Eng.p1.incidence$counts))
## 
## Coefficients:
##        Estimate  Std. Error
## i0 417.30047597 8.265965851
## r   -0.03182554 0.001423115
## 
## -2 log L: -73446.41
```

```r
fit.mle.IoW.p2 <- mle(minusloglik, 
                   start=list(i0=10, r=-0.01), 
                   fixed=list(start.date=start.date.IoW.p2,
                              incidence.dates=IoW.p2.incidence$dates,
                              incidence.counts=IoW.p2.incidence$counts)
                   )
summary(fit.mle.IoW.p2)
```

```
## Maximum likelihood estimation
## 
## Call:
## mle(minuslogl = minusloglik, start = list(i0 = 10, r = -0.01), 
##     fixed = list(start.date = start.date.IoW.p2, incidence.dates = IoW.p2.incidence$dates, 
##         incidence.counts = IoW.p2.incidence$counts))
## 
## Coefficients:
##       Estimate Std. Error
## i0  4.90650182 0.85101721
## r  -0.05320915 0.01030546
## 
## -2 log L: 17.73932
```

```r
fit.mle.UK.p2 <- mle(minusloglik,
                     start=list(i0=10, r=-0.02), 
                     fixed=list(start.date=start.date.UK.p2,
                                incidence.dates=UK.p2.incidence$dates,
                                incidence.counts=UK.p2.incidence$counts))
summary(fit.mle.UK.p2)
```

```
## Maximum likelihood estimation
## 
## Call:
## mle(minuslogl = minusloglik, start = list(i0 = 10, r = -0.02), 
##     fixed = list(start.date = start.date.UK.p2, incidence.dates = UK.p2.incidence$dates, 
##         incidence.counts = UK.p2.incidence$counts))
## 
## Coefficients:
##         Estimate   Std. Error
## i0 1393.11365688 1.447956e+01
## r    -0.01870897 7.094076e-04
## 
## -2 log L: -368702.5
```

