#### Resampling
#
# All resampling techniques take original data, and randomly generate
#  a population of samples from the original data, to get some results.
#  How the population of samples are generated is where most of the
#  differences lie.

############## Permutation test ##############
# The main idea:
#  1. Put the (N1 + N2) results from two groups into a single big bin
#  2. Randomly sample N1 from the big bin into a first group.
#  3. Put the rest into a second group.
#  4. Compare means
#  5. Repeat 1-4 many, many times to get a normal distribution
#  6. Compare the original data to the resulting distribution to get
#     a p-value
# (For the mathematically particular, this should actually be a
#   "combination" test rather than a "permutation" test)

# Example from PLOS One: "Beer consumption increases human attractiveness
#  to malaria mosquitoes" (DOI: 10.1371/journal.pone.0009546). First seen
#  at:
#  http://sas-and-r.blogspot.com/2014/11/example-201413-statistics-doesnt-have.html

beer = c(27, 19, 20, 20, 23, 17, 21, 24, 31, 26, 28, 20, 27, 19, 25, 31, 24, 28, 24, 29, 21, 21, 18, 27, 20)
water = c(21, 19, 13, 22, 15, 22, 15, 22, 20, 12, 24, 24, 21, 19, 18, 16, 23, 20)

ds = data.frame(y = c(beer, water),
                x = c(rep("beer", length(beer)), rep("water", length(water))))

# With mosaic...skip for now
# require(mosaic)
# obsdiff = compareMean(y ~ x, data=ds)
# nulldist = do(999)*compareMean(y ~ shuffle(x), data=ds)
# histogram(~ result, xlab="permutation differences", data=nulldist)
# ladd(panel.abline(v=obsdiff, col="red", lwd=2))

# > obsdiff
# [1] -4.377778
# > tally(~ abs(result) > abs(obsdiff), format="percent", data=nulldist)

# TRUE FALSE
# 0.1  99.9

# Set up the big bin
alldata <- c(beer, water)
labels <- c(rep("beer", length(beer)), rep("water", length(water)))
obsdiff <- mean(alldata[labels=="beer"]) - mean(alldata[labels=="water"])

# Do the permutation tests: Version 1
N <- 100000   # Number of permutations to test
resample_diff <- vector("numeric", N)

# This function only works for the "beer" and "Water" data.
#  It may have to change a bit for general use
resample1 <- function(dat, lab, N) {
  resample_diff <- vector("numeric", N)
  for(i in 1:N) {
    resample_labels <- sample(lab)
    resample_diff[i] <- mean(dat[resample_labels=="beer"]) -
      mean(dat[resample_labels=="water"])
  }
  return(resample_diff)
}

resample_diff <- resample1(alldata, labels, N)
hist(resample_diff)

# Do the permutation tests: Version 2
resamp_means = function(data, labs){
  resample_labels = sample(labs)
  resample_diff = mean(data[resample_labels=="beer"]) -
    mean(data[resample_labels=="water"])
  return(resample_diff)
}
nulldist = replicate(N,resamp_means(alldata,labels))


hist(nulldist, col="cyan")
abline(v = obsdiff, col = "red")
# The p-value is obtained by counting the proportion of statistics
#  (including the actual observed difference) among greater than or equal
#  to the observed statistic:

alldiffs = c(obsdiff,nulldist)
p = sum(abs(alldiffs >= obsdiff)/ N)

############ Bootstrap ##############
#
# Following:
#  http://datascienceplus.com/introduction-to-bootstrap-with-applications-to-mixed-effect-models/
# "[T]he essence of bootstrap: resampling the observed data with
#  replacement and computing the statistic of interest...many times on the
#  resampled data to get a distribution of the statistic of interest.
#  This distribution of the statistic of interest can then be used to
#  compute, for example, confidence intervals."
#
# Bootstraps tend to be robust, except:
#  a) statistic is near the edge of the parameter space. (In such
#      cases, the distribution doesn't converge uniformly, so
#      using the distribution for inferences is problematic)
#  b) bootstrap never converges
#  c) sample size is really small (use resampling instead)
#

# Example

# Generate some data
set.seed(20151101)  # Should always set the seed: Makes results reproducible
height<-rnorm(100,175,6)

# Resampling:
t0<-median(height)
t<-sapply(1:1000,function(x) median(sample(x=height,size=100,replace=TRUE)))
hist(t)
abline(v=t0,col="orange",lwd=3)

# But...the boot package makes it pretty easy to get, for example
#  confidence intervals, too.
library(boot)
b1<-boot(data=height,statistic=function(x,i) median(x[i]),R=1000)
boot.ci(b1)  # boot.ci(b1,conf=0.95) is the default, but other C.I. can
             #  be used, too.
# The original has many other examples, too.

###### Student Course Evaluation Example ######
# Here's an example that is closer to home: Student Course Evaluations
#
# You fill these out, on a scale of 1-7. The instructor gets back
#  all her/his data for each question, along with department and SFJC
#  means, standard deviations, medians, and N. So, for example, in
#  one of my courses, I got the following data for one question:
#               N     Mean    SD    Median
# Instructor    11    6.91    0.30    7
# Department    44    6.59    0.66    7
# College     9279    6.36    1.15    7
#
# I happen to also know that I got 10 @ "7"s and 1 @ "6". (I chose
#  this question because is was the one where I thought I had about
#  the best chance of finding significance.) FWIW, the Department
#  results must be 30 @ "7"s, 10 @ "6"s, and 4 @ "5"s.
# One would be tempted to do an ANOVA or a t-test, because these can
#  be calculated using only the information given, but (a) these both
#  come back not significant, and besides, the distributions are
#  not the correct ones (e.g., normal for an ANOVA).
# Hence, we wrote a C++ program to find all the possible ways that
#  the distribution could be arranged that would give the actual
#  mean, standard deviation, and median. While at it, we also got
#  which of those was the "closest" (in some sense) to the vector
#  distribution of the instructor (i.e.: c(0,0,0,0,0,1,10)).
# Now, we could try something like a "Mann-Whitney" or other non-
#  parametric (rank) test. But, there are too many ties for that to
#  work.
# Hence, the only thing to do: Resample! So here goes:
# We have a distribution, and we need the full thing
inst.dist <- c(0,0,0,0,0,1,10)
# Just some data:
dept.dist <- c(0,0,0,0,4,10,30)
# From the C++ application
## Closest solution:
## 337	1	0	1	0	3909	5031
## Potential solutions: 8,790,293,020
## One close solution
## 0 0 547 0 1786 179 6767
sjfc.dist <- (0,0,547,0,1786,179,6767)
## Another close solution
## 337	1	0	1	0	3909	5031
sjfc.dist <- c(337, 1, 0, 1, 0, 3909, 5031)
ans <- seq(1:7)

# Create the populations from the frequencies
inst.tot <- rep(ans, inst.dist)
dept.tot <- rep(ans, dept.dist)
sjfc.tot <- rep(ans, sjfc.dist)

mean(sjfc.tot)
sd(sjfc.tot)
t0 <- mean(inst.tot)

# Draw 100000 random samples of 11 numbers
num.rep <- 100000
t<-sapply(1:num.rep,
          function(x) mean(sample(x=sjfc.tot,size=11,replace=TRUE)))
hist(t)
abline(v=t0,col="orange",lwd=3)
# p-value:
pv <- length(which(t >= t0)) / num.rep
# Hence, not significant against the Department: p-values ~ 0.08
# Significant against the College, though: p-values ~ 0.012


############## Timing ###############
# Testing the relative timing of some of these:

# Compare!
N <- 10000000
system.time(resample1(alldata, labels, N))
## With N =  1000000: user =  22.81
## With N = 10000000: user = 224.54
system.time(replicate(N,resamp_means(alldata,labels)))
## With N =  1000000: user =  27.96
## With N = 10000000: user = 312.63
# ~20% different, and growing.


############# Jackknife ###########
# This is sort of a cheap bootstrap. In essence, it removes k ("order")
#  observations, and re-computes the quantities of interest to create
#  an estimate of the variance. It does this for each possible set
#  of k observations. (In other words, this is a "leave-k-out"
#  resampling.) Often, k = 1.
#
# It appears that, due to the historical timing, jackknife was
#  quickly replaced by bootstrapping.
#
# The code, unless one uses the "boot" or the "SPECIES" package,
#  is a bit clunky, so we'll pass on it unless you have need of it.
# Here are a couple of references that used it, though. They fall
#  into the time between when jackknife was developed and when
#  bootstrapping because widely available.
#  Burnham, K. P., and Overton,W. S. (1978), Estimation of the
#   Size of a Closed Population When Capture Probabilities Vary
#   Among Animals, Biometrika, 65, 625-633.
#  Burnham, K. P., and Overton,W. S. (1979), Robust Estimation of
#   Population Size When Capture Probabilities Vary Among Animals,
#   Ecology, 60, 927-936.

############# Cross validation #############
# (Some of you may be familiar with "split halves reliability". This
#  is an extension of that.) Essentially, cross validation splits
#  the data set into two or more parts. The first part is used to
#  compute the statistics of interest. Then, those are applied to
#  the second part to see if the two agree. The number of parts
#  is called the number of "folds", and 5-fold or 10-fold are common.
#  In many cases, the cross validation is less biased than the
#  bootstrap analysis, but it is also a bit trickier. Given the ease
#  of bootstrapping, it is probably a good place to start.

package:blocklength (for time series bootstrapping) - see;
https://alecstashevsky.com/r/blocklength/
https://www.r-bloggers.com/2025/02/blocklength-0-2-0/?utm_source=phpList&utm_medium=email&utm_campaign=R-bloggers-daily&utm_content=HTML



