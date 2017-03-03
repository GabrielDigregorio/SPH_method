clc;
close all;
clear all;
set(0,'defaulttextinterpreter','latex');
sz = 50;
%% Load datas
load('khPerformancePert.dat');
load('NperformancePert.dat');
load('Nperformance.dat');
load('khPerformance.dat');

%% Compute proportionality constant
k1 = log(Nperformance(:,2)) - log(Nperformance(:,1).^2);
k1 = exp(mean(k1));
k2 = log(NperformancePert(:,3)) - log(NperformancePert(:,1).*log(Nperformance(:,1)));
k2 = exp(mean(k2));

%% N influence
fig1 = figure(1);
hold on;
grid on;
box('on')
set(gca, 'xscale','log')
set(gca, 'yscale','log')
xlabel('Number of particules $N$ [-]','Fontsize', 12)
ylabel('Time [s]','Fontsize', 12)
set(gca,'xlim',[50 1e4],'Fontsize', 12)

scatter(Nperformance(:,1),Nperformance(:,2),sz,'o','r')
scatter(Nperformance(:,1),Nperformance(:,3),sz,'o','b')

scatter(NperformancePert(:,1),NperformancePert(:,2),sz,'x','r')
scatter(NperformancePert(:,1),NperformancePert(:,3),sz,'x','b')

plot(Nperformance(:,1),k1*Nperformance(:,1).^2,'g')
plot(Nperformance(:,1),k2*Nperformance(:,1).*log(Nperformance(:,1)),'k')

l1 = legend('All-Pair search','Linked-list search','All-Pair search; variable spacing',...
    'Linked-list search; variable spacing','$k_1 N^2$', '$k_2 N$ log$(N)$');
set(l1, 'Interpreter', 'Latex','Fontsize', 12)
set(fig1, 'Position', [100, 100, 650, 350]);
%% kh influence
fig2 = figure(2);
hold on;
grid on;
box('on')
set(gca, 'yscale','log')
xlabel('Smoothing length $\kappa h$ [-]')
ylabel('Time [s]','Fontsize', 12)
set(gca,'xlim',[1.2 3.1],'Fontsize', 12)

scatter(khPerformance(:,1),khPerformance(:,2),sz,'o','r')
scatter(khPerformance(:,1),khPerformance(:,3),sz,'o','b')

scatter(khPerformancePert(:,1),khPerformancePert(:,2),sz,'x','r')
scatter(khPerformancePert(:,1),khPerformancePert(:,3),sz,'x','b')

l2 = legend('All-Pair search','Linked-list search','All-Pair search; variable spacing',...
    'Linked-list search; variable spacing','Fontsize', 12);

set(l2, 'Interpreter', 'Latex','Fontsize', 12)

set(fig2, 'Position', [100, 100, 650, 350]);