function [t,prb]=average_profile(time,probability)

if (time(1)-floor(time(1)))<0.5
    [~,ind1]=min(abs(time-(floor(time(1))+0.5)));
    s=1;
else
    [~,ind1]=min(abs(time-(floor(time(1))+1.5)));
    s=2;
end

if (time(end)-floor(time(end)))<0.5
    e=1;
else
    e=2;
end

l=floor(time(end)-time(1));
if s==2
    l=l-1;
end
if e==1
    l=l-1;
end

t=0.5:0.01:1.5;

for i=1:l
    
    [~,ind2]=min(abs(time-time(ind1)-1));
    time_now=time(ind1:ind2);
    Pr_now=probability(ind1:ind2,:);
    for j=1:size(probability,2)
        prb(:,j,i) = interp1(time_now-floor(time_now(1)),Pr_now(:,j),t);
    end
    ind1=ind2;
end
prb=nanmean(prb,3);