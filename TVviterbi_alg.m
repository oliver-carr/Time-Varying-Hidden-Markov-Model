function [max_ind,delta,delta2,prb]=TVviterbi_alg(obs,Mu,Cov,prior,trans)

prior(prior==0)=1e-200;
p=size(trans,1);

%Compute the observation probabilites, p(x|y)
for i=1:p
    pdf_vals(i,:) = mvnpdf(obs,Mu(i,:),Cov(:,:,i));
end

%1. Initialisation
delta=zeros(p,length(obs(:,1)));
gamma=ones(p,length(obs(:,1)));

%Prevent zeros in the pdf distribution 
ind=pdf_vals==0;
pdf_vals(ind)=0.000001;


for i=1:p
    delta(i,1)=log(prior(1,i))+log(pdf_vals(i,1));
    delta2(i,1)=prior(1,i)*pdf_vals(i,1);    
end
scale(1)=sum(delta2(:,1));
delta2(:,1)=delta2(:,1)./scale(1);

%Recursive calculation
for t=2:length(obs)
    for state=1:p
        delta(state,t)=max(delta(:,t-1)+log(trans(:,state,t)))+log(pdf_vals(state,t));
        delta2(state,t)=max(delta2(:,t-1).*trans(:,state,t))*pdf_vals(state,t);
        [~,max_gamm]=max(delta(:,t-1)+log(trans(:,state,t)));
        prb(state,t)=delta2(:,t-1)'*trans(:,state,t);
        gamma(state,t)=max_gamm;
    end
    scale(t)=sum(delta2(:,t));
    delta2(:,t)=delta2(:,t)./scale(t);
end

%Termination
s_MAP=zeros(1,length(obs(:,1)));
[~,s_MAP(end)]=max(delta(:,end));
Prob_max=max(delta(:,end));

%Backtracking
for t=length(obs)-1:-1:1
    s_MAP(t)=gamma(s_MAP(t+1),t);    
end
max_ind=s_MAP;


















