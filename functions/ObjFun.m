function f = ObjFun(x,T,Z,Xi,N)

% Re-order Xi for matrix multiplication
Xi=permute(Xi,[2 3 1]);


% Z=cos([(time(1)-median(diff(time)))-x(1,1,end);time-x(1,1,end)]*2*pi);
% Z=[cos([(time(1)-median(diff(time)));time]*2*pi),sin([(time(1)-median(diff(time)));time]*2*pi)];

% 
% for t=1:T-1
%     tempT(t)=squeeze(sum(sum(Xi(:,:,t).*log(exp(x(:,:,1)+x(:,:,2)*Z(t))./sum(exp(x(:,:,1)+x(:,:,2).*Z(t)),2)),2),1));
% end
% f=1/sum(tempT);

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
