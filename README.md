# Hidden Markov Model Toolbox for Matlab

Matlab implementation of standard hidden Markov models (HMMs) with either discrete observations and continuous emissions, and dependent HMMs which allow the parameters to vary with time. Further details and examples are provided below for both supervised and unsupervised cases with multivariate data:

# Hidden Markov Model (Discrete)

## Supervised case

COMING SOON

## Unsupervised case

COMING SOON

# Continuous Emission Hidden Markov Model

Inputs:
* X - an M by T matrix of observed continuous data. M is the number of observed variables and T is the total time, observations must be evenly sampled but may contain missing data.

* N - the number of hidden states

Outputs:
* S - sequence of hidden states. 1 by T vector of most likely hidden states from the Viterbi algorithm

* Pi - probabilities of initial states. 1 by N vector of probabilities.

* mu - means of the normal distributions. N by M vector of means for each observation in each state.

* cov - covariance matrices of the normal distributions. M by M by N matrix of covariances for each observation in each state.

* A - transition matrix of probabilities. N by N matrix of transitions between each state.

### Example

```matlab
ContinuousEmissionHMM.m
```

#### Example 1 - Estimation of sleep/wake (rest/active) from single total acceleration observations (N=2)

![Total acceleration data and transformed data (log differences)](/images/contData1.png)
Total acceleration data and transformed data (log differences).

![Data Distribution](/images/contDist1.png)
Gaussian distributions for the two hidden states (active and inactive) after running the Baum-Welch algorithm.

![States](/images/contStates1.png)
Viterbi sequence of hidden states with observation data, and the probability of each state at each time point. Inactive state is marked by the blue area and active state by the red area.

#### Example 2 - Estimation of sleep/wake (rest/active) from the three orthogonal acceleration components which make up Example 1 (N=2)

![Acceleration in the three directions](/images/contData2.png)
Acceleration data in each of the three orthogonal directions.

![Probability of states](/images/contStates2.png)
Probability of the inactive state (red) and active state (blue) at each timepoint using multivariate acceleration observations.



### Time Varying Hidden Markov Model
