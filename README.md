# Hidden Markov Model Toolbox for Matlab

Matlab implementation of standard hidden Markov models (HMMs) with either discrete observations and continuous emissions, and dependent HMMs which allow the parameters to vary with time. Further details and examples are provided below for both supervised and unsupervised cases with multivariate data.

### Functions
* Baum-Welch algorithm - function to determine HMM parameters from unsupervised set of observations.

* Viterbi algorithm - function to determine the most likely sequence of hidden states from a set of observations and HMM parameters.

### Parameters
* State Sequences - a set of sequences of observed (non-hidden) states, each sequence may be of variable length.

* Observation Sequences - a set of sequences of observations corresponding, each sequence contains O observations and may be of variable length. Observations may be continuous or discrete.

* N - the number of hidden states.

* Pi - probabilities of initial states.

* ContObs - parameters which define the probability distribution of continuous observations in each state (mean and covariance for normal distribution).

* DiscObs - probabilities of discrete observations in each state.

* A - transition matrix of probabilities. N by N matrix of transitions between each state.

* S - sequence of hidden states. Sequence of most likely hidden states from the Viterbi algorithm for each of the observation sequences for testing.


# Example - Supervised Hidden Markov Model (Discrete)



# Example - Unsupervised Hidden Markov Model (Discrete)



# Example - Supervised Hidden Markov Model (Continuous Emission)

# Example - Unsupervised Hidden Markov Model (Continuous Emission)

[Example code](./ContinuousEmissionHMMexample.m) which attempts to determine active and inactive periods from accelerometry data. First using a single continuous observation, total acceleration, to predict active and inactive periods. Second, three continuous observations, the three orthogonal components of acceleration are used to predict active and inactive periods as well as a third intermediate activity period. 

### Example 1 - Estimation of sleep/wake (rest/active) from single total acceleration observations (N=2)

![Total acceleration data and transformed data (log differences)](/images/contData1.png)
Total acceleration data and transformed data (log differences).

![Data Distribution](/images/contDist1.png)
Gaussian distributions for the two hidden states (active and inactive) after running the Baum-Welch algorithm.

![States](/images/contStates1.png)
Viterbi sequence of hidden states with observation data, and the probability of each state at each time point. Inactive state is marked by the blue area and active state by the red area.

### Example 2 - Estimation of sleep/wake (rest/active) from the three orthogonal acceleration components which make up Example 1 (N=2)

![Acceleration in the three directions](/images/contData2.png)
Acceleration data in each of the three orthogonal directions.

![Probability of states](/images/contStates2.png)
Probability of the inactive state (red) and active state (blue) at each timepoint using multivariate acceleration observations.

![Probability of three states](/images/cont3States2.png)
Three state HMM with three orthogonal acceleration observations. Probability of the inactive state (yellow), partially active state (red) and active state (blue) at each timepoint using multivariate acceleration observations.



# Time Varying Hidden Markov Model
