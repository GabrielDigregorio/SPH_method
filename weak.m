%% weak scaling
close all 
clear all
clc
%% first test
% timetot=[807.669 955.742 1103.71 1056.07 938.013 941.655 930.642 969.594 971.426 979.737 1222.97 2138.95];
% number_process=[1 2 3 4 5 6 8 10 15 20 40 80];
% timeintegration=[801.309 910.974 969.637 919.532 849.12 851.405 828.83 846.347 838.281 834.768 867.746 1451.31];
% timewriting=[3.04529 41.1788 130.46 132.113 84.1559 86.2661 97.7158 118.785 127.192 138.839 346.991 674.3];
% timeinitwritting=[0.312322 0.538456 0.678553 1.58242 0.866747 1.03127 1.33148 1.59272 2.31718 3.05732 5.28579 10.7528];
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% 
% plot(number_process,timetot,number_process,timeintegration,number_process,timewriting,number_process,timeinitwritting,'LineWidth',4)
% hold on 
% h = legend('Total time','Time of $timeintegration$','Time of wirrting and exhanging','Time init(exchan+init+writt0)');
% set(h,'Interpreter','latex','Location','NorthWest')
% ylim([-10 2200])
% title(' ','Interpreter','latex','FontSize',40)
% set(gca,'fontsize',40)
% hold on 
% grid on
% xlabel(' Number of CPU','Interpreter','latex')
% ylabel('Computation time','Interpreter','latex')
% set(gcf,'color','w')

%% second test 
% timetot=[758.581 869.425 1057.89 921.616 1019.16  926.897  1256.01  1240.06  963.037 983.954 1043.66 1187.26];
% number_process=[1 2 3 4 5 6 8 10 15 20 40 80];
% timeintegration=[752.685 824.437 916.928 836.037 899.544 837.784 945.541 940.561 829.355 832.409 823.547 821.806];
% timewritingloop=[2.95215 6.50413 9.78822 13.0716 16.7512 19.4418 26.641 34.4502 48.6694 65.0537 130.441 263.762];
% timeMPIloop=[0.00297284 35.5629 128.081 69.285 98.7568 65.8811 279.915 260.311 79.5028 80.3371 80.6954 88.1127];
% 
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% 
% plot(number_process,timetot,number_process,timeintegration,number_process,timewritingloop,number_process,timeMPIloop,'LineWidth',4)
% hold on 
% h = legend('Total time','Time of $timeintegration$','Time of wirrting in loop','Time MPI');
% set(h,'Interpreter','latex','Location','NorthWest')
% ylim([-10 1500])
% title(' ','Interpreter','latex','FontSize',40)
% set(gca,'fontsize',40)
% hold on 
% grid on
% xlabel(' Number of CPU','Interpreter','latex')
% ylabel('Computation time','Interpreter','latex')
% set(gcf,'color','w') 




%% third test 
number_process=[1 2 3 4 5 6 8 10 15 20 40 80];
timetot=[733.128 877.04 1114.07 903.489 967.853 908.655 978.16 959.766 1109.28 1012.44 1075.25 1185.93];
timeintegration=[ 727.581 827.726 975.98 821.792 873.154 820.841 851.629 846.305 939.357 861.765 837.479 821.345];
timewritingloop=[2.92897 6.53317 10.2176 13.0286 16.4106 19.4079 25.8061 32.5071 51.2258 65.4506 130.745 263.106];
timeMPIloop=[0.00236058 39.8416 124.285 64.9126  74.2188 64.7517 95.3952 76.0271 113.471 79.0725 99.1609 87.9291];
timeMPIintegration=[0.00161242 10.4628 89.1693 37.3003 41.8205 36.4368 64.8696 47.2154 81.5876 47.424 66.3704 49.3419];
timeinit=[0.286288 0.496931 0.61546 0.744577 0.919904 1.0794 1.28017 1.58121 2.33027 3.08445 5.30438 10.3767];
timewritinginit=[0.122498 0.257484 0.410681 0.471249 0.6003 0.707717 0.93552 1.18244 1.83807 2.38491 4.20675 8.56473];

figure('units','normalized','outerposition',[0 0 1 1])
hold on 
plot(number_process,timetot,'-o','color',[0 0.5 0],'LineWidth',2.5)
plot(number_process,timeMPIintegration,'-o','color',[0 0 0.5],'LineWidth',2.5)

plot(number_process,timeintegration-timeMPIintegration,'-o','color',[0.5 0 0],'LineWidth',2.5)%
plot(number_process,timeMPIloop,'-o','color',[0.5 0.5 0],'LineWidth',2.5)
plot(number_process,timewritingloop+timewritinginit,'-o','color',[0.5 0.25 0],'LineWidth',2.5)
%plot(number_process,timeinit-timewritinginit,'-o','color',[0.5 0.25 0.25],'LineWidth',2.5)

h = legend('Total','Communication in \textit{timeintegration}','Time of numerical integration ','Communication each time step','Writing time');
set(h,'Interpreter','latex','Location','NorthWest')
ylim([-10 1250])
title(' ','Interpreter','latex','FontSize',40)
set(gca,'fontsize',40)
hold on 
grid on
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 20]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
xlabel(' Number of processor $p$[-]','Interpreter','latex','FontSize',28)
ylabel('Computation time[s]','Interpreter','latex','FontSize',28)
set(gcf,'color','w')


%% in /100

figure('units','normalized','outerposition',[0 0 1 1])
hold on 
plot(number_process,timetot./timetot*100,'-o','color',[0 0.5 0],'LineWidth',2.5)
plot(number_process,timeMPIintegration./timetot*100,'-o','color',[0 0 0.5],'LineWidth',2.5)

plot(number_process,(timeintegration-timeMPIintegration)./timetot*100,'-o','color',[0.5 0 0],'LineWidth',2.5)%
plot(number_process,timeMPIloop./timetot*100,'-o','color',[0.5 0.5 0],'LineWidth',2.5)
plot(number_process,(timewritingloop+timewritinginit)./timetot*100,'-o','color',[0.5 0.25 0],'LineWidth',2.5)


h = legend('Total','Communication in $timeintegration$','Time of numerical integration ','Communication each time step','Writing time');
set(h,'Interpreter','latex','Location','NorthWest')
ylim([-5 110])
title(' ','Interpreter','latex','FontSize',40)
set(gca,'fontsize',40)
hold on 
grid on
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 20]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
xlabel(' Number of processor $p$[-]','Interpreter','latex','FontSize',28)
ylabel('Computation time[\%]','Interpreter','latex','FontSize',28)
set(gcf,'color','w')


%% sequential and paralle
figure('units','normalized','outerposition',[0 0 1 1])
hold on 
plot(number_process,timetot./timetot*100,'-o','color',[0 0.5 0],'LineWidth',2.5)
plot(number_process,timeMPIintegration./timetot*100,'-o','color',[0 0 0.5],'LineWidth',2.5)

plot(number_process,(timeintegration-timeMPIintegration)./timetot*100,'-o','color',[0.5 0 0],'LineWidth',2.5)%
plot(number_process,timeMPIloop./timetot*100,'-o','color',[0.5 0.5 0],'LineWidth',2.5)
plot(number_process,(timewritingloop+timewritinginit)./timetot*100,'-o','color',[0.5 0.25 0],'LineWidth',2.5)


h = legend('Total','Communication in $timeintegration$','Time of numerical integration ','Communication each time step','Writing time');
set(h,'Interpreter','latex','Location','NorthWest')
ylim([-5 110])
title(' ','Interpreter','latex','FontSize',40)
set(gca,'fontsize',40)
hold on 
grid on
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 20]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
xlabel(' Number of processor $p$[-]','Interpreter','latex','FontSize',28)
ylabel('Computation time[\%]','Interpreter','latex','FontSize',28)
set(gcf,'color','w')