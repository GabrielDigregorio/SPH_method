%% Plot the position and the speed of the wavefront of the Dam_break over the time.

%% Change Experiment_name depending on the name of your dam_break experiment in the folder 'Results'. 
Experiment_name = 'dam';
NAME = strcat('../build/Results/',Experiment_name);

[time,t,X,t1_exp,X1_exp,t2_exp,X2_exp,t3_exp,X3_exp,v] = front(NAME);


figure(1)
    hold on
    plot(t(1:22),X(1:22))
    plot(t1_exp,X1_exp,'o')
    plot(t2_exp,X2_exp,'s')
    plot(t3_exp,X3_exp,'rx')
    set(gca,'xlim',[0,3.5])
    set(gca,'ylim',[1 4])
    grid on
    box on
    ylabel('$Y/L$','Interpreter','Latex')
    xlabel('$t(2g/L)^{1/2}$','Interpreter','Latex')
    set(gca,'FontSize',25)
    h=legend('Numerical results','Experimental results');
    set(h,'Interpreter','Latex')
    
    figure(2)
    hold on
    plot(time,v)
    grid on
    box on
    ylabel('$u_y$ [m/s]','Interpreter','Latex')
    xlabel('Time [s]','Interpreter','Latex')
    set(gca,'FontSize',25)
  