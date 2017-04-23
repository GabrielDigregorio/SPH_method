%% Strong scaling analysis
% Case - Dam break - Dam3
% Free particles: 32175
% BC particles: 24546

processors = [1 2 3 4 5 6 8 10 12 15 17 20 25 30];
contProc = linspace(1, 30, 4);

time = [48.98 28.93 18.71 13.15 10.82 8.82 6.62 5.33 4.7 4.01 3.45 2.82 2.82 2.08];


figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 20]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on
plot(processors, time(1)./time, '-o', 'linewidth', 2.5);
plot(contProc, contProc, 'k', 'linewidth', 1);


leg = legend('Code speedup','Ideal speedup $S=p$',...
    'Location','southeast');
set(leg,'Interpreter','latex')

xlabel('Number of processors $p$ [-]','Interpreter','latex','FontSize',28);
ylabel('Speedup $S=T_1/T_p$ [-]','Interpreter','latex','FontSize',28);
%xlim([-0.1,0.6]);
%ylim([-1,1.5]);

