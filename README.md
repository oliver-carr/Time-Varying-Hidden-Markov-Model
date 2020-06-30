# Hidden Markov Model Toolbox for Matlab

Matlab implementation of standard hidden Markov models (HMMs) with continuous emissions, and dependent HMMs which allow the parameters to vary with time. The time varying hidden Markov model has been developed from the work of [Huang et al.](https://royalsocietypublishing.org/doi/full/10.1098/rsif.2017.0885) and applied in [Monitoring Depression in Bipolar Disorder using Circadian Measures from Smartphone Accelerometers](https://scholar.google.co.uk/citations?user=dpFEilMAAAAJ&hl=en).

When using this code please cite:
[1] O. Carr, F. Andreotti, K. E. A. Saunders, N. Palmius, G. M. Goodwin, M. De Vos, “Monitoring Depression in Bipolar Disorder using Circadian Measures from Smartphone Accelerometers”.
[2] Q. Huang, D. Cohen, S. Komarzynski, X. M. Li, P. Innominato, F. Lévi, and B. Finkenstädt, “Hidden Markov models for monitoring circadian rhythmicity in telemetric activity data,” J. R. Soc. Interface, vol. 15, no. 139, 2018.

### Functions
* Baum-Welch algorithm - function to determine HMM parameters from unsupervised set of observations.
* Time Varying Transition Probability Baum-Welch - Baum-Welch algorithm with time varying transition probabilities.

* Viterbi algorithm - function to determine the most likely sequence of hidden states from a set of observations and HMM parameters.
* Time Varying Viterbi Algorithm - function to determine the most likely sequence of hidden states from a set of observations and time varying HMM parameters.
 
### Parameters
* State Sequences - a set of sequences of observed (non-hidden) states, each sequence may be of variable length.

* Observation Sequences - a set of sequences of observations corresponding, each sequence contains O observations and may be of variable length. Observations may be continuous or discrete.

* N - the number of hidden states.

* Pi - probabilities of initial states.

* A - transition matrix of probabilities. N by N matrix of transitions between each state, or N by N by T of time varying transitions between each state.

* S - sequence of hidden states. Sequence of most likely hidden states from the Viterbi algorithm for each of the observation sequences for testing.


# Example Outputs

[Example code](./ContinuousEmissionHMMexample.m) which attempts to determine active and inactive periods from accelerometry data. First using a single continuous observation, total acceleration, to predict active and inactive periods. Second, three continuous observations, the three orthogonal components of acceleration are used to predict active and inactive periods as well as a third intermediate activity period. Finally, a single continuous observation is used with a time varying covariate which allows the transition probabilities to vary with time.

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

### Example 3 - Estimation of sleep/wake (rest/active) from single total acceleration observations (N=2) with a time varying covariate

Two time varying series are used as covariates. A sine wave and cosine wave with 24 hour period to ensure no repeating coviarate values for each 24 hour period.

![Total acceleration data and transformed data with viterbi (log differences)](/images/Time Varying Data.png)
Total acceleration data and transformed data (log differences) with viterbi sequence.

![Probability of states and average weekly profile](/images/Time Varying States.png)
Probability of the inactive and active states at each time point, with the average weekly profile below.

