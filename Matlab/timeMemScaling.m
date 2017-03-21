N=10.^(1:1:7);
timeNaive=[0,0,0,0,0,0,0];
timeList=[0,0,0,0,0,0,0];
memNaive=[0,0,0,0,0,0,0];
memList=[0,0,0,0,0,0,0];

%% Time plot
    figure;
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[0.0 0.0 60 60*3/4]);
    set(gcf,'PaperPosition',[0.0 0.0 60 60*3/4]);
    grid on;
    box on;
    set(gca,'FontSize',22);
    xlabel('$N$ (log)','FontSize',22,'Interpreter','latex');
    ylabel('Time [s]','FontSize',22,'Interpreter','latex');
    hold on;
    %plot(N,timeNaive,'b.-');
    %plot(N,timeList,'r.-');
    semilogx(N,timeNaive,'b.-');
    semilogx(N,timeList,'r.-');
    hold off;
    %axis([0 4 -1 1]);
    NumTicks=5;
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    legendInfo={'All-pair';'Linked List'};
    legend(legendInfo,'Interpreter','latex','Location','Best');
    
%% Mem plot
    figure;
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[0.0 0.0 60 60*3/4]);
    set(gcf,'PaperPosition',[0.0 0.0 60 60*3/4]);
    grid on;
    box on;
    set(gca,'FontSize',22);
    xlabel('$N$ (log)','FontSize',22,'Interpreter','latex');
    ylabel('Memory []','FontSize',22,'Interpreter','latex');
    hold on;
    plot(N,memNaive,'b.-');
    plot(N,memList,'r.-');
    hold off;
    %axis([0 4 -1 1]);
    NumTicks=5;
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    legendInfo={'All-pair';'Linked List'};
    legend(legendInfo,'Interpreter','latex','Location','Best');