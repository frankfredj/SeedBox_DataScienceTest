#What is the aproximate probability distribution between the test group and the control group

Both variables exhibit the same probabilistic structure, namely:

⋅⋅* They have a certain probability of making transaction(s) after visiting the cancelation page
⋅⋅* Given they do perform transaction(s), they will spawn a vector {Revenue, ReBill, ChargeBack, ReFund} 

with mean vector µ, and variance-covariance matrix Σ. 

Sampling from a numerical pdf f(Revenue, ReBill, ChargeBack, ReFund) is quite cumbersome in nature, as 

all variables are dependent. 

From a quick glimpse at the numerical pdf ([Link to plots](https://imgur.com/a/U3A0qvF)), variables are 

not normaly distributed at all. The amount of clients that have made transaction(s) isn't that big, too: 

only 1079 for the control group, and 1635 for the test one. Common assumtions in regards to transactions 

is that they are Poisson distributed. This would give rise to a compound Poisson distribution for 

revenues, and our sample size is too small to give relevant credibility to any estimator {µ, σ} derived 

under a normality assumption. (A common credibility standard is n/λ = 1082.41 (1 + σ^2 / µ^2) for 

Compound Poisson processes.)

I have my own bias in regards to assuming normaly distributed µ's, because I do not know how big "n" 

needs to be. This is why I suggest using Monte Carlo simulations to transform the problem into a 

strictly binomial one. 

In order to do so, I propose the following method:

⋅⋅* Diagonalize Σ with the transformation matrix Ω, thus finding indepedent linear combinations of 

{Revenue, ReBill, ChargeBack + ReFund}. (I've added the negative transactions together because of their 

sparsity). 
⋅⋅* Estimate the probability density function of the linear combinations {PC1, PC2, PC3}
⋅⋅* Build a function that simulates {PC1, PC2, PC3} in accordance to their respective pdf's
⋅⋅* Simulate {PC1, PC2, PC3}, and recover {Revenue, ReBill, ChargeBack + ReFund} with the inverse of Ω

This circumvent the issues that arise from covariances while preserving the fundamental structure of the 

sample's {µ, Σ}

Hence, for example, a random variable taken from the {Test} distribution would be (0,0,0) with 

probability 0.8897877, and equal to Ω^(1) (PC1, PC2, PC3)^T with probability 0.11021234. 

Note that this extend the possible values of {Revenue, ReBill, ChargeBack + ReFund} from {R, N, N} to 

{R, R, R}. This, however, isn't a problem if we use our simulation function to calculate likelihoods.


All of this might seem a bit too fancy / complex, but it has it's advantages: it will give rise to a 

method that is much more intuitive than abstract t-tests slapped with normality assumptions.

Note: The code file also includes a PCA graph and a logistic regression on transactions. This was done 

in order to highlight the failure to distinguish between a {Test} and a {Control} variable based on 

predictive modeling with respect to transactions. All in all, trying to model some function f

(transaction) with output (Test, Control) wouldn't be wise.


#Is a user that must call-in to cancel more likely to generate at least 1 addition REBILL?

We'll simmulate 2 vectors of random variables, namely:

1. A difference between {REBILL | test} and {REBILL | control} labelled R(test, control)
2. A difference between {REBILL | control} and {REBILL | control} labelled R(control, control)

We will then compute respective p's, where "p" is the proportion of R >= 1.

Since p's are techically derived from a sum of Bernouli random variables, their summation exhibit a 

Binomial distribution with parameters (p, n). Such a distribution converges rapidly towards a Normal 

distribution N(p, p(1-p)/n). We've effectively cleared out the problem of having a sufficiently large 

"n" where the underlying distribution of the random variable is unknown in order to apply a normality 

assumption, namely for t-tests.

The idea behind all of this is to simply check what are the average odds of observing a difference of 

more than 1 REBILL between two people who can opt-out online. If implementing a different opt-out system 

really does affect the odds described above, then our newly estimated "p" will be an outlier.

Using the simulation tools described in the previous section with n = 100 000, we have gathered the 

following statistics:

⋅⋅* p given {test, control} was equal to 0.07309 and had variance 6.774785e-07 given n
⋅⋅* p given {control, control} was equal to 0.01920 and had variance 1.883136e-07 given n

Our hypothesis is that: p given {test, control} - p given {control, control} = 0

We can safely label the distribution of this random variable as normal with mean 0 and s.d. 6.774785e-07 

+ 1.883136e-07 under our hypothesis.

The resulting p-value is 1, which indicates that someone amongst the test group is indeed likely to 

generate at least 1 more rebill. (far-right outlier)



#Is a user that must call-in to cancel more likely to generate more revenues?

Working in the exact same manner as earlier, we just swap the REBILL variable for the total REVENUE one. 

Also, instead of calculating the proportion of differences greater than or equal to 1, we'll simply 

calculate the proportion of difference greater than 0.

Using n = 100 000 again, we have gathered the following statistics:

⋅⋅* p given {test, control} was equal to 0.10545 and had variance 9.43303e-07 given n
⋅⋅* p given {control, control} was equal to 0.02414 and had variance 2.355726e-07 given n

The resulting p-value is 1, which indicates that someone amongst the test group is indeed likely to 

generate more revenue. (far-right outlier)


#Is a user that must call-in more likely to produce a higher chargeback rate(CHARGEBACKs/REBILLs)?

For this metric, we'll use {Revenue, Chargeback / ReBill, Refund} as our base variables. Also, we will 

use the following rules:

⋅⋅* Chargeback / ReBill = 1 if ReBill = 0, but Chargeback > 0
⋅⋅* Chargeback / ReBill = 0 if ReBill = 0 and Chargeback = 0

Working in same Binomial settings with n = 100 000, we have gathered the following statistics:

⋅⋅* p given {test, control} was equal to 0.10545 and had variance 9.43303e-07 given n
⋅⋅* p given {control, control} was equal to 0.02414 and had variance 2.355726e-07 given n

The resulting p-value is 0, which indicates that someone amongst the test group not likely to have a 

higher chargeback rate. (far-left outlier)


#Conclusion

Based on Monte Carlo simulations of Bernouli Random Variables respecting the fundamental structure {µ, 

Σ} of our sample, we can safelly declare that forcing people to opt-out via a phone-in system is likely 

to:

⋅⋅* Generate more ReBill
⋅⋅* Increase revenues
⋅⋅* Decrease ChargeBack rates

