function [Mu,Cov,A,Pi,x]=TVTP_HMM(X,time,N,cyc)

rng(3)
lik=0;
LL=[];

p=length(X(1,:));
T=length(X(:,1));

if nargin<5   tol=0.00001; end
if nargin<4   cyc=100; end
if nargin<3   N=2; end

% Set the initial covariate
Z=[cos([(time(1)-median(diff(time)));time]*2*pi),sin([(time(1)-median(diff(time)));time]*2*pi)];

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
A=repmat(A,1,1,T);

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
            alpha(t,i)=sum(alpha(t-1,:).*A(:,i,t)')*B(t,i);
        end
        scale(t)=sum(alpha(t,:));
        alpha(t,:)=alpha(t,:)./scale(t);
    end
    
    Scale=Scale+log(scale');
    
    % Beta recursion with same scaling
    beta(T,:)=ones(1,N)/scale(T);
    for t=T-1:-1:1
        for i=1:N
            beta(t,i)=sum(A(i,:,t).*B(t+1,:).*beta(t+1,:))/scale(t);
        end
    end
    
    % E (Expectation step) Estimate the state occupation probabilities
    for t=1:T-1
        for i=1:N
            for j=1:N
                Xi(t,i,j)=alpha(t,i)*A(i,j,t)*B(t+1,j)*beta(t+1,j);
            end
        end
        a=Xi(t,:,:);
        Xi(t,:,:)=a./sum(a(:));
    end
    Gamma=(sum(Xi,3));
    
    
    % M (Maximization step) Re-estimate the HMM parameters
    Pi=Gamma(1,:);
    
    % CALCULATE A(i,j,t)
    x10=rand(N);    
    x0(:,:,1)=x10;
    for cv=1:size(Z,2)
        x0(:,:,cv+1)=rand(N);
    end
    x0(:,:,cv+2)=zeros(size(N,N));
    
    fun=@(x)ObjFun(x,T,Z,Xi,N);
    
    lb = -Inf*ones(N,N,size(Z,2)+1);
    ub = Inf*ones(N,N,size(Z,2)+1);
    for i=1:N
        lb(i,i,:)=0;
        ub(i,i,:)=0;
    end 
    lbEnd=zeros(N,N);
    ubEnd=zeros(N,N);
    lbEnd(1,1) = -0.5;
    ubEnd(1,1) = 0.5;
    
    lb(:,:,end+1)=lbEnd;
    ub(:,:,end+1)=ubEnd;
    
    Acon = [];
    bcon = [];
    Aeq = [];
    beq = [];
    nonlcon=[];
    
%     options=optimoptions('fmincon','Display','none','ConstraintTolerance',1e-5,'Algorithm','active-set','MaxIterations',200,'OptimalityTolerance',1e-5,'StepTolerance',1e-5);
    
    tic
    x = fmincon(fun,x0,Acon,bcon,Aeq,beq,lb,ub);
    toc
    
    Z=sin([(time(1)-median(diff(time)))-x(1,1,end);time-x(1,1,end)]*2*pi);
    
    for t=1:T
        for i=1:N
            for j=1:N
                A(i,j,t)=exp(x(i,j,1)+x(i,j,2)*Z(t))/sum(exp(x(i,:,1)+x(i,:,2).*Z(t)));
%                 A(i,j,t)=exp(x(i,j,1)+squeeze(x(i,j,2:end))'*Z(t,:)')/sum(exp(x(i,:,1)+(squeeze(x(i,:,2:end))*Z(t,:)')'));
            end
        end
    end
    
    
    
    for i=1:N
        Mu(i,:)=sum(X(1:end-1,:).*Gamma(:,i))./sum(Gamma(:,i));
        
        d=X(1:end-1,:)-Mu(i,:);
        Cov(:,:,i)=((d.*(Gamma(:,i)*ones(1,p)))'*d)./sum(Gamma(:,i));
    end
    %     Cov=sqrt(Cov);
    
    oldlik=lik;
    lik=sum(Scale);
    LL=[LL lik];
    fprintf('cycle %i log likelihood = %f ',cycle,lik);
    if (cycle<=2)
        likbase=lik;
    elseif (lik<oldlik)
        fprintf('violation');
    elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase)||~isfinite(lik))
        fprintf('\n');
        break;
    end
    fprintf('\n');
    
    
end
