# What is the aproximate probability distribution between the test group and the control group

Both variables exhibit the same probabilistic structure, namely:

* They have a certain probability of making transaction(s) after visiting the cancelation page
* Given they do perform transaction(s), they will spawn a vector {Revenue, ReBill, ChargeBack, ReFund} with mean vector µ, and variance-covariance matrix Σ. 

From a quick glimpse at the numerical pdf ([Link to plots](https://imgur.com/a/U3A0qvF)), variables are not normaly distributed at all. The amount of clients that have made transaction(s) isn't that big, too: only 1079 for the control group, and 1635 for the test one. Common assumtions in regards to transactions is that they are Poisson distributed. This would give rise to a compound Poisson distribution for revenues, and our sample size is too small to give relevant credibility to any estimator {µ, σ} derived under a normality assumption. (A common credibility standard is n/λ = 1082.41 (1 + σ^2 / µ^2) for Compound Poisson processes.)

I have my own bias in regards to assuming normaly distributed µ's, because I do not know how big "n" needs to be. This is why I suggest using a variant of Monte Carlo simulations instead. In order to do so, I propose the following method:


* Diagonalize Σ with the transformation matrix Ω, thus finding independent linear combinations of {Revenue, ReBill, ChargeBack + ReFund}. (I've added the negative transactions together because of their sparsity). 
* Estimate the probability density function of the linear combinations {PC1, PC2, PC3}
* Build a function that simulates {PC1, PC2, PC3} in accordance to their respective pdf's
* Simulate {PC1, PC2, PC3}, and recover {Revenue, ReBill, ChargeBack + ReFund} with the inverse of Ω

This circumvent the issues that arise from covariances while preserving the fundamental structure of the sample's {µ, Σ}

Hence, for example, a random variable taken from the {Test} distribution would be (0,0,0) with probability 0.8897877, and equal to Ω^(1) (PC1, PC2, PC3)^T derived in accordance to the test sample's {µ, Σ} with probability 0.11021234. (0.11021234 is the proportion of the test group that went on to complete transactions after being labelled.)

(Note that this extend the possible values of {Revenue, ReBill, ChargeBack + ReFund} from {R, N, N} to {R, R, R}. This, however, isn't a problem if we use our simulation function to calculate likelihoods.)

All of this might seem a bit too fancy / complex, but it has it's advantages: it will give rise to a method that is much more intuitive than abstract t-tests slapped with normality assumptions.

Note: The code file also includes a PCA graph and a logistic regression on transactions. This was done in order to highlight the failure to distinguish between a {Test} and a {Control} variable based on predictive modeling with respect to transactions. All in all, trying to model some function f(transaction) with output (Test, Control) wouldn't be wise.




# Is a user that must call-in to cancel more likely to generate at least 1 addition REBILL?

We'll simmulate 2 vectors of random variables, namely:

1. A difference between {REBILL | test} and {REBILL | control} labelled R(test, control)
2. A difference between {REBILL | control} and {REBILL | control} labelled R(control, control)

