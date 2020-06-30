function f = ObjFun(x,T,Z,Xi,N)

% Re-order Xi for matrix multiplication
Xi=permute(Xi,[2 3 1]);

c=0;
for t=1:T-1
    for j=1:N
        for k=1:N
            c=c+1;
            temp(c)=Xi(j,k,t)*log(exp(x(j,k,1)+squeeze(x(j,k,2:end))'*Z(t,:)')/sum(exp(x(j,:,1)+(squeeze(x(j,:,2:end))*Z(t,:)')')));
        end
    end
end
f=1/sum(temp);
