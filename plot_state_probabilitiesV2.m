function plot_state_probabilitiesV2(prb,time)

figure

h=area(time,prb','LineStyle','none');
h(1).FaceColor = [0 0 1];
h(2).FaceColor = [1 0.4 0.4];
% h(3).FaceColor = [1 0 0];
datetick
axis tight
xlabel('Date/Time','fontsize',14)
ylabel('Probability of each state','fontsize',14)
