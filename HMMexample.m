slashchar = char('/'*isunix + '\'*(~isunix));
mainpath = (strrep(which(mfilename),[mfilename '.m'],''));
addpath(genpath(mainpath)) % add subfunctions folder to path

% Load the example file
[Acc_Time,Total_Acc,Acc,HR_Time,HR]=load_example_data();

% Define the data to be used in the HMM
Example=1; %Select Example to run

if Example==1
    %Example 1 - Transformed total acceleration to find 3 hidden states of
    %activity (inactive, mildy active, highly active)
    time=Acc_Time(2:end);
    data=log(abs(diff(Total_Acc)));
    S=2; %Number of hidden states
    
elseif Example==2
    %Example 2 - Three-dimensional activity to find 3 hidden states of
    %activity (inactive, mildy active, highly active)
    time=Acc_Time;
    data=Acc;
    S=3; %Number of hidden states
    
elseif Example==3
    %Example 3 - Heart rate to find two hidden states
    time=HR_Time;
    data=HR;
    S=3; %Number of hidden states
end

% Remove any -Inf values from log(0) and replce with minimum + noise
data(data==-Inf)=NaN;
nans=find(isnan(data));
min_d=min(data);
for nn=1:length(nans)
    data(nans(nn))=min_d+(nanstd(data)/2)*randn(1,1);
end

% Plot the data (For Example 1)
figure('units','normalized','outerposition',[0 0 1 1])
subplot(5,1,1)
plot(Acc_Time,Total_Acc,'k.')
datetick
axis tight
title('HMM Example','fontsize',15)
ylabel('Raw Data','fontsize',14)
subplot(5,1,2)
for pl=1:size(data,2)
    hold on
    plot(time,data(:,pl),'.')
end
datetick
axis tight
ylabel('Data','fontsize',14)
% Plot the distribution of data
subplot(5,1,3)
for pl=1:size(data,2)
    hold on
    histogram(data(:,pl),50,'Normalization','pdf')
end
ylabel('Distribution','fontsize',14)

% Calculate mean, covariance, tranistion matrix and prior using the forward
% backward algorithm (Baum-Welch)
[Mu,Cov,A,Pi]=forward_backward(data,S);

% Plot the pdfs on the histogram (only for 1D input)
if size(data,2)==1
    [~,indMu]=sort(Mu(:,1));
    Mu=Mu(indMu,:);
    Cov=Cov(:,:,indMu);
    A=A(indMu,:);
    Pi=Pi(:,indMu);
    
    for i=1:S
        pd{i} = makedist('Normal',Mu(i),Cov(:,:,i));
        hold on
        plot(linspace(min(data),max(data),10000),pdf(pd{i},linspace(min(data),max(data),10000)))
        dist=pdf(pd{i},linspace(min(data),max(data),1000));
        pctval = prctile(dist,95);
    end
end

% Find the hidden states at each observation and the probabilities of each
% state
[max_ind,delta,prb]=viterbi_alg(data,Mu,Cov,Pi,A);

% Plot the state transitions
plot_state_probabilities(prb,time,data,max_ind,S)











