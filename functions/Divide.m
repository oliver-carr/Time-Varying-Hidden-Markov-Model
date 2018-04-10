function [idx]=Divide(X,K,type) 
%X is data, K is number of groups to divide into

if strcmp(type,'sort')
    [~,I] = sort(X);
elseif strcmp(type,'rand')
    I = randperm(length(X));
else
    I = randperm(length(X));
end

T=length(X);

for i=1:K
    idx{i}=I(floor((i-1)*T/K +1):floor(i*T/K));
end