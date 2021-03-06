# What is the aproximate probability distribution between the test group and the control group

Both variables exhibit the same probabilistic structure, namely:

* They have a certain probability of making transaction(s) after visiting the cancelation page
* Given they do perform transaction(s), they will spawn a vector {Revenue, ReBill, ChargeBack, ReFund} with mean vector µ, and variance-covariance matrix Σ. 

Given that a user did perform a transactions, the numerical probability density functions look like this (linear interpolation was used to fit both distribution over the same domain, rather than using bar plots):

![](https://i.imgur.com/IEMxlnt.png)

And the ReBill, ChargeBack, Refund look like:

![](https://i.imgur.com/uQq6h3j.png)

To showcase the differences in Revenue, ReBill, ChargeBack and Refund amongst the two groups, we can refer the the Principal Component plot:

![](https://i.imgur.com/1DnZOoO.png)

From a quick glimpse at the graphs above, variables are not normaly distributed, and both groups are not easily separable. The amount of clients that have made transaction(s) isn't that big, too: only 1079 for the control group, and 1635 for the test one. Common assumtions in regards to transactions is that they are Poisson distributed. This would give rise to a compound Poisson distribution for revenues, and our sample size with regards to transactions is too small to give relevant credibility to any estimator {µ, σ} derived under a normality assumption. (A common credibility standard is n/λ = 1082.41 (1 + σ^2 / µ^2) for Compound Poisson processes.)

I have my own bias in regards to assuming normaly distributed µ's, because I do not know how big "n" needs to be. This is why I suggest using a variant of Monte Carlo simulations rather than classic t-tests on means. In order to do so, I propose the following method:


* Diagonalize Σ with the transformation matrix Ω, thus finding independent linear combinations of {Revenue, ReBill, ChargeBack + ReFund}. (I've added the negative transactions together because of their sparsity). 
* Estimate the probability density function of the linear combinations {PC1, PC2, PC3}
* Build a function that simulates {PC1, PC2, PC3} in accordance to their respective pdf's
* Simulate {PC1, PC2, PC3}, and recover {Revenue, ReBill, ChargeBack + ReFund} with the inverse of Ω

This circumvents the issues that arise from sampling given non-zero covariances while preserving the fundamental structure of the sample's {µ, Σ}. (In my code file, I give an example of 100 000 random values being simulated from the {test | transaction} sample. The resuting mean vector and variance-covariance matrix is virtualy identical to the sample's mean vector and variance-covariance matrix. The output is shown below:)

![](https://i.imgur.com/4D5rY2b.png)

Hence, for example, a random variable taken from the {Test} distribution would be (0,0,0) with probability 0.8897877, and equal to Ω^(-1) (PC1, PC2, PC3)^T derived in accordance to the test sample's {µ, Σ} with probability 0.11021234. (0.11021234 is the proportion of the test group that went on to complete transactions after being labelled.)

So, for any question that is binomial in nature (i.e.: are test subjects likely to produce at least one more ReBill), we could derive the following comprehensible test:

1. Simulate differences between 2 random variables drawn from the {Control} sample. 
2. Simulate differences between 2 random variables drawn from the {Test} and {Control} sample.
3. Compute the respective likelyhoods, that is, the proportion of differences that are greater than a certain treshold. (For the ReBill case, this treshold would be 1.)
4. Check if the {Test - Control} proportion is an outlier with regards to the {Control - Control} distribution

Assuming normality of the {Control - Control} is reasonable since it is a Binomial variable with parameters (n,p), which converges rapidly to a normal distribution.

We do, however, need to draw a variable from the samples where a transaction has occured in order to obtain a non-zero difference. So our Normal distribution with regards to {Control} will have estimated mean p, and variance p(1-p)/1079.



# Is a user that must call-in to cancel more likely to generate at least 1 addition REBILL?

We'll simmulate 2 vectors of random variables, namely:

1. A difference between {REBILL | test} and {REBILL | control} labelled R(test, control)
2. A difference between {REBILL | control} and {REBILL | control} labelled R(control, control)

We will then compute respective p's, where "p" is the proportion of R >= 1.

Using the simulation tools described in the previous section with n = 100 000, we have gathered the following statistics:

* p given {test, control} was equal to 0.07309 
* p given {control, control} was equal to 0.01920 

Our hypothesis is that: p given {test, control} - p given {control, control} = 0

We can now perform a z-test with variance =  0.01920(1-0.01920)/1079 and mean = 0.01920 to check how likely it would be to observe p = 0.07309. 

The above test results give a p-value of 1, indicating that users amongst the test groups are likely to generate at least one more ReBill than their control counterparts. (Note that 7% of the consumers are expected to generate one more ReBill. This is a small increase, but a statisticaly significant one nonetheless)



# Is a user that must call-in to cancel more likely to generate more revenues?

Working in the exact same manner as earlier, we just swap the REBILL variable for the total REVENUE one. 

Also, instead of calculating the proportion of differences greater than or equal to 1, we'll simply calculate the proportion of differences greater than 0.

Using n = 100 000 again, we have gathered the following statistics:

* p given {test, control} was equal to 0.10545 
* p given {control, control} was equal to 0.02414 

The resulting p-value is 1, which indicates that someone amongst the test group is indeed likely to generate more revenue (since it is a far-right outlier). 



# Is a user that must call-in more likely to produce a higher chargeback rate(CHARGEBACKs/REBILLs)?

For this metric, we'll use {Revenue, Chargeback / ReBill, Refund} as our base variables. Also, we will use the following rules:

* Chargeback / ReBill = 1 if ReBill = 0, but Chargeback > 0
* Chargeback / ReBill = 0 if ReBill = 0 and Chargeback = 0

Working in same Binomial settings with n = 100 000, we have gathered the following statistics:

* p given {test, control} was equal to 0.02632
* p given {control, control} was equal to 0.02429 

The resulting p-value is 0.667545, which is inconclusive. Perhaps the outcome would have been different if we just checked the difference for Chargeback rather than Chargeback / ReBill. But under these settings, we reject the hypothesis that a user who must call-in is more likely to produce a higher chargeback rate.



# Conclusion

Based on Monte Carlo simulations of Bernouli Random Variables respecting the fundamental structure {µ, Σ} of our sample, we can safely declare that forcing people to opt-out via a phone-in system is likely to:

* Generate more ReBill (Around 7% of the consumers will generate more ReBill)
* Increase revenues (Around 10% of the consumers will generate more revenues)



