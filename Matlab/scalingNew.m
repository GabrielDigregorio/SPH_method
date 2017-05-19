%% Strong scaling analysis

% GOOD (MPI and OpenMP)
% Case - Dam break - Dam3
% Free particles: 32175
% BC particles: 24546

% BAD (MPI only)
% Case - Dam on wet bed - Dam4
% Free particles: 37800
% BC particles: 24546

%% MPI
processors = [1 2 3 4 5 8 10 12 15 17];
contProc = linspace(1, 17, 4);

timeGood = [89.5; 
    46.94;
    34.25;
    25.6;
    22.3;
    12.5;
    10;
    9.9;
    7.75;
    8.33];

timeBad = [79.2; 
    59.3;
    60;
    49.7;
    52.8;
    28.5;
    23.5;
    21.5;
    20.1;
    19.8];


figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 20]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on
plot(processors, timeGood(1)./timeGood, '-o', 'linewidth', 2.5);
plot(processors, timeBad(1)./timeBad, '-o','color', 'm', 'linewidth', 2.5);
plot(contProc, contProc, 'k', 'linewidth', 1);


leg = legend('Speedup for the ideal case', 'Speedup for the non-ideal case', '$S=p$',...
    'Location','north');
set(leg,'Interpreter','latex')

xlabel('Number of processors $p$ [-]','Interpreter','latex','FontSize',28);
ylabel('Speedup $S=T_1/T_p$ [-]','Interpreter','latex','FontSize',28);
xlim([0,17]);
%ylim([-1,1.5]);


%% OpenMP
processorsOMP = [1 2 3 4 5 6 8 10 13 16];
contProcOMP = linspace(1, 16, 4);
% With 2 nodes
time2 = [47.5; 
    25.35;
    17.1;
    12.5;
    10.8;
    9.25;
    7.3;
    7.3;
    5.5;
    5];
% With 5 nodes
time5 = [20.12; 
    12.5;
    7.25;
    6;
    5.8;
    4.3;
    4;
    4;
    3.4;
    3.1];


figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 20]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on
plot(processorsOMP, time2(1)./time2, '-o', 'linewidth', 2.5);
plot(processorsOMP, time5(1)./time5, '-o','color', 'm', 'linewidth', 2.5);
plot(contProc, contProc, 'k', 'linewidth', 1);


leg = legend('Speedup with 2 nodes', 'Speedup with 5 nodes', '$S=CPU_n$',...
    'Location','north');
set(leg,'Interpreter','latex')

xlabel('Number of CPU per node $CPU_n$ [-]','Interpreter','latex','FontSize',28);
ylabel('Speedup $S=T_1/T_{CPU_n}$ [-]','Interpreter','latex','FontSize',28);
xlim([0,16]);
%ylim([-1,1.5]);



