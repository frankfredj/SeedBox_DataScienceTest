# What is the aproximate probability distribution between the test group and the control group

Both variables exhibit the same probabilistic structure, namely:

* They have a certain probability of making transaction(s) after visiting the cancelation page
* Given they do perform transaction(s), they will spawn a vector {Revenue, ReBill, ChargeBack, ReFund} with mean vector µ, and variance-covariance matrix Σ. 

From a quick glimpse at the numerical pdf ([Link to plots](https://imgur.com/a/U3A0qvF)), variables are not normaly distributed at all. The amount of clients that have made transaction(s) isn't that big, too: only 1079 for the control group, and 1635 for the test one. Common assumtions in regards to transactions is that they are Poisson distributed. This would give rise to a compound Poisson distribution for revenues, and our sample size is too small to give relevant credibility to any estimator {µ, σ} derived under a normality assumption. (A common credibility standard is n/λ = 1082.41 (1 + σ^2 / µ^2) for Compound Poisson processes.)
