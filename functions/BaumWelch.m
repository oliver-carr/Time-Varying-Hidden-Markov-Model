function [Mu,Cov,trans,Pi]=BaumWelch(X,N,cyc,tol)
% Gaussian Observation Hidden Markov Model
%
%  Input:
%       - X:                Observations (MxT)
%       - N:                Number of states
%       (optional inputs)
%           - cyc:          Number of iterations of the forward-backward algorithm   
%           - tol:          Tolerance to stop iterations
%           
%   Output:
%       - Mu:                Means of the normal distributions for each
%                            observation in each state (MxN)
%       - Cov:               Covariance of the normal distributions for eachs
%                            observation in each state (MxMxN)
%       - Pi:                Initial state probabilities (1xN)
%       - trans:             Transition matrix (NxN)
% 
% --
% Released under the GNU General Public License
%
% Copyright (C) 2019  Oliver Carr
% University of Oxford, Insitute of Biomedical Engineering, CIBIM Lab - Oxford 2017
% fernando.andreotti@eng.ox.ac.uk
%
% 
% For more information visit: https://github.com/fernandoandreotti/cinc-challenge2017
% 
% Referencing this work
%
% Last updated : April 2019
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


lik=0;
LL=[];

p=length(X(1,:));
T=length(X(:,1));

if nargin<2
    error('Not enough input arguments')
elseif nargin==2
    cyc=100;
    tol=0.00001;
elseif nargin==3
    tol=0.00001;
end


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
trans=rand(N);
trans=trans./sum(trans')';

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
            alpha(t,i)=sum(alpha(t-1,:).*trans(:,i)')*B(t,i);
        end
        scale(t)=sum(alpha(t,:));
        alpha(t,:)=alpha(t,:)./scale(t);        
    end
    
    Scale=Scale+log(scale');
    
    % Beta recursion with same scaling
    beta(T,:)=ones(1,N)/scale(T);
    for t=T-1:-1:1
        for i=1:N
            beta(t,i)=sum(trans(i,:).*B(t+1,:).*beta(t+1,:))/scale(t);
        end
    end
    
    % E (Expectation step) Estimate the state occupation probabilities
    for t=1:T-1
        for i=1:N
            for j=1:N
                Xi(t,i,j)=alpha(t,i)*trans(i,j)*B(t+1,j)*beta(t+1,j);
            end
        end
        a=Xi(t,:,:);
        Xi(t,:,:)=a./sum(a(:));
    end    
    Gamma=(sum(Xi,3));
    
    
    % M (Maximization step) Re-estimate the HMM parameters
    Pi=Gamma(1,:);
    
    for i=1:N
        trans(i,:)=squeeze(sum(squeeze(Xi(:,i,:)),1))/squeeze(sum(Gamma(:,i)));
    end
    
    
    for i=1:N
        Mu(i,:)=sum(X(1:end-1,:).*Gamma(:,i))./sum(Gamma(:,i));       
        
        d=X(1:end-1,:)-Mu(i,:);        
        Cov(:,:,i)=((d.*(Gamma(:,i)*ones(1,p)))'*d)./sum(Gamma(:,i));
    end
end



















