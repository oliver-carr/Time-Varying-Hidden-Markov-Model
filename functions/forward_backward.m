% Gaussian Observation Hidden Markov Model

% Implementing EM and Viterbi algorithms for Hidden Markov Model
% in linear memory
%
% A set of states S={S_1,S_2,...,S_N}, with q_t being the state visited at
% time t
%
% A set of PDFs B={b_1(o),...,b_N(o)}, describing the emission
% probabilities b_j(o_t)=p(o_t|q_t=S_j) for 1<=j<=N where o_t is the
% observations at time point t from the sequence O={o_1,o_2,...,o_T}}
%
% The state transition probability matrix A={a_ij} for 1<=i,j<=N where
% a_ij=p(q_t+1=S_j|q_t=S_i)
%
% The initial state distribution vector Pi={pi_1,...,pi_N}
%
% lamba=(Pi,A,B) are a set of parameters completely specifying the HMM
%
% Forward procedure
% alpha_t(i)=p(o_1,...,o_t|q_t=S_i,lambda)
% Initially alpha_1(i)=pi_i b_i(o_1) for 1<=i<=N
%
% alpha_t(j)=[sum(i=1 to N) alpha_t-1(i) a_i,j] b_j(o_t) for t=2,3,...,T
% and 1<=j<=N
%
% p(O|lambda)=sum(i=1 to N) alpha_T(i) is the sequence likelihood
%
% Backward procedure
% beta_t(i)=p(o_t+1,...,o_T|q_t=S_i,lambda)
% Initially beta_T(i)=1 for 1<=i<=N
%
% beta_t(i)=sum(j=1 to N) a_i,j b_j(o_t+1) beta_t+1(j) for t=T-1,...,1 and
% 1<=i<=N
%
% p(O|lambda)=sum(i=1 to N) pi_i b_i(o_1) beta_1(i)
%
% Xi_t(i,j) is defined as the probability of being in state i at time t,
% and state j at time t + 1, given the model and the observation sequence
%
% Xi_t(i,j)=p(q_t=S_i,q_t+1=S_j|O,lambda)
%
% Xi_t(i,j)=alpha_t(i) a_i b_j(o_t+1) beta_t+1(j) / p(O|lambda)
%
% Gamma_t(i) is defined as the probability of being in state i at time t,
% given the observation sequence and the model
%
% Gamma_t(i)=p(q_t=S_i|O,lambda)
%
% Gamma_t(i)=alpha_t(i) beta_t(i) / sum(i=1 to N){alpha_t(i) beta_t(i)}
%
% Gamma_t(i)=sum(j=1 to N) {Xi_t(i,j)}
%
% X - T x p data matrix
% N - number of states (default 2)
% cyc - maximum number of cycles of Baum-Welch (default 100)
%
% Mu - mean vectors
% Cov - output covariance matrix (full, tied across states)
% A - state transition matrix
% Pi - priors
%
% Iterates until cyc steps of Baum-Welch

function [Mu,Cov,A,Pi]=forward_backward(X,N,cyc)

rng(3)
lik=0;
LL=[];

p=length(X(1,:));
T=length(X(:,1));

if nargin<4   tol=0.00001; end
if nargin<3   cyc=100; end
if nargin<2   N=2; end

% Find the initial means and covariance matrices for each of the states
% Split the observations into evenly size states from smallest to largest
[idx]=Divide(X,N,'sort');

Cov=zeros(p,p,N);
for i=1:N
    Cov(:,:,i)=cov(X(idx{i},:));
    Mu(i,:)=mean(X(idx{i},:));
end

% Initialise Priors
Pi=rand(1,N);
Pi=Pi/sum(Pi);

% Initialise Transition matrix
A=rand(N);
A=A./sum(A')';

for cycle=1:cyc
    
    %%%% FORWARD-BACKWARD
    Scale=zeros(T,1);
    
    % Find the probabilty of each observation being from each state
    for i=1:N
        B(:,i) = mvnpdf(X,Mu(i,:),Cov(:,:,i));
    end
    B(B==0)=1e-200;
    
    % Initial alpha step
    alpha(1,:)=Pi.*B(1,:);
    scale(1)=sum(alpha(1,:));
    alpha(1,:)=alpha(1,:)./scale(1);
    
    % Alpha recursion with scaling to avoid underflow
    for t=2:T
        for i=1:N
            alpha(t,i)=sum(alpha(t-1,:).*A(:,i)')*B(t,i);
        end
        scale(t)=sum(alpha(t,:));
        alpha(t,:)=alpha(t,:)./scale(t);        
    end
    
    Scale=Scale+log(scale');
    
    % Beta recursion with same scaling
    beta(T,:)=ones(1,N)/scale(T);
    for t=T-1:-1:1
        for i=1:N
            beta(t,i)=sum(A(i,:).*B(t+1,:).*beta(t+1,:))/scale(t);
        end
    end
    
    % E (Expectation step) Estimate the state occupation probabilities
    for t=1:T-1
        for i=1:N
            for j=1:N
                Xi(t,i,j)=alpha(t,i)*A(i,j)*B(t+1,j)*beta(t+1,j);
            end
        end
        a=Xi(t,:,:);
        Xi(t,:,:)=a./sum(a(:));
    end    
    Gamma=(sum(Xi,3));
    
    
    % M (Maximization step) Re-estimate the HMM parameters
    Pi=Gamma(1,:);
    
    for i=1:N
        A(i,:)=squeeze(sum(squeeze(Xi(:,i,:)),1))/squeeze(sum(Gamma(:,i)));
    end
    
    
    for i=1:N
        Mu(i,:)=sum(X(1:end-1,:).*Gamma(:,i))./sum(Gamma(:,i));       
        
        d=X(1:end-1,:)-Mu(i,:);        
        Cov(:,:,i)=((d.*(Gamma(:,i)*ones(1,p)))'*d)./sum(Gamma(:,i));
    end
end



















