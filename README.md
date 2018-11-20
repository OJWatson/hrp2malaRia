# hrp2malaRia

Individual based malaria transmission model for the study of *P. falciparum* *pfhrp2* gene deletions. The model is an extension of the Imperial College transmission model used as part of the [Global Technical Strategy 2030 costing exercise](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(15)00423-5/fulltext), and was used to characterise the drivers of *pfhrp2* deletion selection in the following publication:

 > Watson OJ, Slater HC, Verity R, et al. Modelling the drivers of the spread of Plasmodium falciparum hrp2 gene deletions in sub-Saharan Africa. Elife 2017; 6: e25008. (https://elifesciences.org/articles/25008)

Installing *hrp2malaRia* devel
-------------
To install the development version from github (the package *devtools is required*):

```r
library(devtools)
install_github("OJWatson/hrp2malaRia")
```

Once installed, the package can be loaded using:

```r
library(hrp2malaRia)
```

Asking a question
------------------
- for bug reports, feature requests, contributions, use github's [issue system](https://github.com/OJWatson/hrp2_malaRia/issues)

## Demonstration

We can demonstrate this by simulating 2 populations of 2000 individuals. Each 
population starts with 10% frequency of hrp2 deletions. One population only
uses RDTs, whereas the other uses 70% RDTs, as well as not adhering to a
negative RDT result and still treating patients. Both populations have an 
annual EIR of 1. 

```r

# population 1
r <- hrp2_Simulation(N=2000, years=20, rdt.det = 0, EIR=1/365)

# population 2
r2 <- hrp2_Simulation(N=2000,years=20,rdt.det = .3,microscopy.use = .3, EIR = 1/365)

# and let's plot these - all veriables that begin with S. are series variables collected over time. 

# let's plot the first population
plot(r$S.Times,r$S.N.Dels, xlab = "Time (days)", ylab = "hrp2 Deletion Frequency", ylim=c(0,1), col="red", type="l", lwd=2)

# and add the second population
lines(r2$S.Times,r2$S.N.Dels,lwd=2)

# and add our legend
legend(1500, .15, legend=c("Only RDT use", "70% RDT use, and 30% nonadherence to RDT negative result"),
col=c("red", "black"), lty=1, cex=0.8,
box.lty=2, box.lwd=2)
```
![](tools/demo.png)

## Input Variables to Change

| Variable | Description                                                                   |
|----------|-------------------------------------------------------------------------------|
| *years* | Number of years simulation is run for. 20 years will be usually long enough to see hrp2 selection at low EIRs. |
| *N* | Population size. The model is stochastic and so a larger N will enable clearer patterns to be found in model outputs, however, they will take longer to run. Population of 2000 should be a good balance of speed and clarity for demonstration purposes. For the analysis in the eLife modelling paper N was set to 100,000 |
| *strains.0* | Starting ratio of strains in the population. This number must be an integer. The default is set to 10, which will mean that on initialisation there will be roughly 10 times as many wild type strains as there are hrp2 deleted. |
| *EIR* | Annual EIR, provided as a fraction, i.e. the default EIR of 20/365 represents an EIR of 20. |
| *ft* | Frequency of treatment being sought upon developing clinical disease. Must be a number between 0 and 1, with the default being 0.4, i.e. 40% of cases seek treatment. |
| *rdt.det* | The probability that an individual only infected with hrp2-deleted strains still produces a positive RDT. The default is 1, i.e. no selective pressure. For the eLife paper we used a value of 0.25 that represented the approximate 25% chance that a positive RDT would still be produced due to hrp3 protein epitope cross reactivity |
| r*dt.nonadherence* | The probability that the result of an RDT is ignored, i.e. a negative RDT due to hrp2 deletions will still be treated. Default = 0, and in the eLife paper we used a value of 10% based on data collected by the WHO. |
| *microscopy.use* | The proportion of cases that are tested by microscopy rather than RDT. Default = 0, and in the eLife paper we tested the impact of 30% of cases being diagnosed by microscopy, representing the ~70% of cases diagnosed being diagnosed by RDTs in 2014. |
| *fitness* | Comparative fitness of hrp2 deleted parasites. Default = 1, which mean the wild type is equally as likely to be passed on to mosquitoes as an hrp2 deleted parasite. 0.8 results in the hrp2 deleted strain being 80% likely to be passed on, whereas a value of 1.2 would mean it is 20% more likely. |

## Output Variables of Interest

| Variable | Description                                                                   |
|----------|-------------------------------------------------------------------------------|
| *S.Times* | The timing of when the simulated population is sampled, i.e. the time at which all other series variables (variables that begin S.) were collected |
| *S.Prev.All* | The prevalence of malaria within the whole population |
| *S.Prev.05* | The prevalence of malaria within under 5s |
| *S.N.Dels* | The proportion of all parasites that are hrp2 deleted |
| *S.N.Dels.05* | The proportion of all parasites within under 5s that are hrp2 deleted |
| *S.Incidence* | The daily clinical incidence of malaria. Multiply this by 365 to have the expected number of clinical cases of malaria a person within the whole population is expected to have in a year |
| *S.Incidence.05* | The daily clinical incidence of malaria. Multiply this by 365 to have the expected number of clinical cases of malaria a child under the age of 5 is expected to have in a year |
| *S.Prev.Mono.D* | The proportion of infected individuals who are only infected with hrp2 deleted parasites |
| *S.Prev.Mono.D.05* | The proportion of infected individuals under the age of 5 who are only infected with hrp2 deleted parasites |
