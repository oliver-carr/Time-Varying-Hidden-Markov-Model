# Hidden Markov Model Toolbox for Matlab

Matlab implementation of standard hidden Markov models (HMMs) with either discrete observations and continuous emissions, and dependent HMMs which allow the parameters to vary with time. Further details and examples are provided below for both supervised and unsupervised cases with multivariate data:

### Hidden Markov Model (Discrete)

* Supervised case

COMING SOON

* Unsupervised case

COMING SOON

### Continuous Emission Hidden Markov Model

Inputs:
* X - an M by T matrix of observed continuous data. M is the number of observed variables and T is the total time, observations must be evenly sampled but may contain missing data.

* N - the number of hidden states

Outputs:
* S - sequence of hidden states. 1 by T vector of most likely hidden states from the Viterbi algorithm

* Pi - probabilities of initial states. 1 by N vector of probabilities.

* mu - means of the normal distributions. N by M vector of means for each observation in each state.

* cov - covariance matrices of the normal distributions. M by M by N matrix of covariances for each observation in each state.

* A - transition matrix of probabilities. N by N matrix of transitions between each state.

Example

```matlab
ContinuousEmissionHMM.m
```

* Example 1 - Estimation of sleep/wake (rest/active) from single total acceleration observations (N=2)






### Time Varying Hidden Markov Model
