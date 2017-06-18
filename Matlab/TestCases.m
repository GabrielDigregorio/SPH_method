%**************************************************************************
% Test Cases:
%                 - Free falling cube under gravity
%                 - Stationary tank
%                 - Crashed Cube
%**************************************************************************
function sucess = TestCases(nameExperiment, n, path)

%close all;
set(groot,'defaultLineLineWidth',2)

% Check path (optional arument)
if nargin < 2
   path = '';
end

%addpath(genpath('../build/Results/'))

    % Number of files to read
    dirName = dir([nameExperiment, '/*.txt']); %list all directory with.txt
    nstep = length(dir([nameExperiment, '/*.txt']))-1 %[-]

    % Import data at t=0
    filename=strcat(nameExperiment,'/',dirName(1).name);
    InitExperiment = importdata(filename);
    
    % time step
    Str = char(InitExperiment.textdata(9));Key = 'Step Time (k) :';
    Index = strfind(Str, Key);
    timeStep = sscanf(Str(Index(1) + length(Key):end), '%g', 1); %[s]
    
    % write step
    Str = char(InitExperiment.textdata(10));Key = 'Write interval : ';
    Index = strfind(Str, Key);
    timeWrite = sscanf(Str(Index(1) + length(Key):end), '%g', 1); %[s]
    
    % Simulation time
    Str = char(InitExperiment.textdata(11));Key = 'Simulation Time (T) : ';
    Index = strfind(Str, Key);
    timeSimu = sscanf(Str(Index(1) + length(Key):end), '%g', 1); %[s]

    % Current Simulation time
    Str = char(InitExperiment.textdata(12));Key = 'Current Time Simulation : ';
    Index = strfind(Str, Key);
    currTimeSimu = sscanf(Str(Index(1) + length(Key):end), '%g', 1); %[s]
    
    % Domain and nbr particules
    Str1 = char(InitExperiment.textdata(13));
    Key1 = 'Domain (lower l) : ';
    Str2 = char(InitExperiment.textdata(14));
    Key2 = 'Domain (upper u) : ';
    Str3 = char(InitExperiment.textdata(15));
    Key3 = 'Number of Particules (nFree/nMoving/nFixed) : ';             
    Index1 = strfind(Str1, Key1); Index2 = strfind(Str2, Key2); Index3 = strfind(Str3, Key3);
    
    % Domain (l and u)
    l     = sscanf(Str1(Index1(1) + length(Key1):end), '%g %g %g', 3);
    u     = sscanf(Str2(Index2(1) + length(Key2):end), '%g %g %g', 3);
    
    % Number of particules (nFree, nMoving and nFixed)
    limit = sscanf(Str3(Index3(1) + length(Key3):end), '%g %g %g', 3); 
    
    % figure(100)
    %     plot(InitExperiment.data(:,3), InitExperiment.data(:,8))

    for i=1 : nstep

        % Open File nbr i
        %disp(dirName(i+1).name);
        filename=strcat(nameExperiment,'/',dirName(i).name);
        Experiment = importdata(filename); % Import Data
        
        % Compute time of the experiment
        time(i) = (i-1)*timeWrite; % Time 
        
        % Get time consumption
        Str = char(Experiment.textdata(6)); Key = 'CPU Time : ';
        Index = strfind(Str, Key);
        CPU_Time(i) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
        
        % Get memory consumption
        Str1 = char(Experiment.textdata(7)); Key1 = 'Memory Usage :';
        Str2 = char(Experiment.textdata(8)); Key2 = 'Peak :';
        Index1 = strfind(Str1, Key1); Index2 = strfind(Str2, Key2);
        Memory(i)      = sscanf(Str1(Index1(1) + length(Key1):end), '%g', 1);
        Memory_Peak(i) = sscanf(Str2(Index2(1) + length(Key2):end), '%g', 1);
    
    end
    
    
    % Put file in a new directory
    %     if(input('MoveFile to an other folder (1 for yes, 0 for no): ' ))
    %     mkdir ../build/Results/FreeFallingCube
    %     movefile ../build/Results/FreeFallingCube_*.txt ../build/Results/FreeFallingCube
    %     movefile ../build/Results/FreeFallingCube_*.vtk ../build/Results/FreeFallingCube
    %     end



switch n


%% Free Falling cube
%  ************************************************************************
case 1

    
    % Parameters
    g = 9.81;
    z0_center = 10; %[m]
    
    for i=1 : nstep 
        % Open File nbr i
        %disp(dirName(i+1).name);
        filename=strcat(nameExperiment,'/',dirName(i).name);
        Experiment = importdata(filename); % Import Data
        
        % Compute the error to analytical solution 
        Analytic(i) = mean(InitExperiment.data(:,3)) - g*(time(i)^2)/2; % MRUA
        MeanZ_experiment(i) = mean(abs(Experiment.data(:,3)));
        error(i) = abs((MeanZ_experiment(i)-Analytic(i))./Analytic(i) *100); % error [%]
        %XY_move(i) = sqrt(mean(  (InitExperiment.data(:,1) - Experiment.data(:,1)).^2  ...
                             % +  (InitExperiment.data(:,2) - Experiment.data(:,2)).^2   ));
    end
 
    
    % Plot
    figure(1)
    hold on
    plot(time, Analytic, time, MeanZ_experiment);
        title('')
        xlabel('Time [s]')
        ylabel('Z')
        text =  strcat('Simulation (time step = ', num2str(timeStep), ')');
        legend('Analytic', text)
        grid
        %print('FreeFallingCube_error', '-depsc')
    hold off  
    

    % Plot
    figure(2)
    hold on
    plot(time, error, '*');
        title('')
        xlabel('Time [s]')
        ylabel('Error in %')
        legend('')
        grid
        %print('FreeFallingCube_error', '-depsc')
    hold off  
    
    figure(3)
    hold on
    plot(time, Memory*1.25e-7, time, Memory_Peak*1.25e-7);
        title('Memory Consumption')
        xlabel('Time [s]')
        ylabel('Memory Consumption [MB]')
        legend('RSS Memory', 'RSS Memory Peak')
        grid
        %print(strcat(path,'FreeFallingCube_Memory'), '-depsc')
    hold off 
    
    figure(4)
    hold on
    plot(time, CPU_Time);
        title('Time Consumption')
        xlabel('Time Experiment[s]')
        ylabel('CPU Time [s]')
        grid
        %print(strcat(path,'FreeFallingCube_Memory'), '-depsc')
    hold off 

% Save data
DATA.name = dirName(1).name;
DATA.path = nameExperiment;
DATA.nbrFiles = nstep;
DATA.timeStep = timeStep;
DATA.timeWrite = timeWrite;
DATA.timeSimu = timeSimu;
DATA.currTimeSimu = currTimeSimu;
DATA.CPUtime = CPU_Time;
DATA.memory = Memory;
DATA.memoryPeak = Memory_Peak;
DATA.time = time;
DATA.analytic = Analytic;
DATA.meanZexperiment = MeanZ_experiment;
DATA.error = error;

save(strcat(nameExperiment,strcat('/',dirName(1).name(1:end-13))), 'DATA')


%% Swimming Pool (Hydrostatic)
%  ************************************************************************
case 2

    % Parameters
    nbrWindows = 6;
    % Cube:


    for i=1 : nstep 
        
        % Open File nbr i
        %disp(dirName(i+1).name);
        filename=strcat(nameExperiment,'/',dirName(i).name);
        Experiment = importdata(filename); % Import Data
        
        % Zmin and Zmax
        Zmax = max(Experiment.data(1:limit(1),3));
        Zmin = min(Experiment.data(1:limit(1),3));
        
        % Height of the swimming pool
        Height(i) = Zmax - Zmin;
        
        % Sort by height
        [Experiment.data(1:limit(1),3), SortIndex] = sort(Experiment.data(1:limit(1),3));
        for j=[1 2 4 5 6 7 8 9]
            Experiment.data(1:limit(1),j) = Experiment.data(SortIndex,j);
        end
        

    
        % Compute the hydrostatic pressure
        for j=1:nbrWindows
            height_min = (j-1)*(Height(i)/nbrWindows);%(nbrWindows-(j-1))*(Height/nbrWindows);
            height_max = (j)*(Height(i)/nbrWindows);%(nbrWindows-(j))*(Height/nbrWindows);
            WindowsDown = min(find(Experiment.data(1:limit(1),3) >= height_min));
            WindowsUp = max(find(Experiment.data(1:limit(1),3) <= height_max));
            
            H(i,j) = height_min + (height_max - height_min)/2;%mean(Experiment.data(WindowsDown:WindowsUp,3));
            mean_Density(i,j)= mean(Experiment.data(WindowsDown:WindowsUp,7));
            std_Density(i,j)= std(Experiment.data(WindowsDown:WindowsUp,7));
            mean_Hydrostatic(i,j)= mean(Experiment.data(WindowsDown:WindowsUp,8));
            std_Hydrostatic(i,j)= std(Experiment.data(WindowsDown:WindowsUp,8));
        end
        
    end


% Delta Log
[pks,locs] = findpeaks(mean_Hydrostatic(1:end/2,3),time(1:end/2));
Delta_log = mean(log(pks(2:end-1)./pks(3:end)));
epsilon_log = sqrt(Delta_log^2/(Delta_log^2+4*pi^2))

% natural Frequency of the system
Deltat = mean(locs(2:end)-locs(1:end-1));
freq = 1/Deltat

figure(1)
hold on
    for i=[1 floor(length(H(:,1))/4) floor(length(H(:,1))/2) 3*floor(length(H(:,1))/4) length(H(:,1))-1]
        errorbar(H(i,:), mean_Density(i,:),std_Density(i,:),'LineWidth', 2)
    end   
    %axis([0 1 2 3])
    title('Density')
    xlabel('Height [m]')
    ylabel('Density [kg/m^3]')
    legend('t = 0 [s]', 't = end/4 [s]', 't = end/2 [s]', 't = 3*end/4 [s]', 't = end [s]')
    grid
    %print(strcat(path,'FreeFallingCube_Memory'), '-depsc')
    hold off   

figure(2)
hold on
    plot(H(1,:), mean_Hydrostatic(1,:), '*','LineWidth', 2)
    %i=[floor(length(H(:,1))/4) floor(length(H(:,1))/2) 3*floor(length(H(:,1))/4) length(H(:,1))-1]
    for i=[floor(length(H(:,1))/4)  length(H(:,1))]
        errorbar(H(i,:), mean_Hydrostatic(i,:), std_Hydrostatic(i,:), '*','LineWidth', 2)
    end
    plot(H(1,:), 1000*9.81*(Height(end)+0.025 -H(1,:)), '-k')
    %plot(H(end,:), 1279.7*9.81*(Height(end) -H(end,:)), '-.k')
    %axis([0 1 2 3])
    set(gca,'FontSize',22);
    xlabel('Height [m]','FontSize',22,'Interpreter','latex');
    ylabel('Hydrostatic Pressure [Pa]','FontSize',22,'Interpreter','latex');
    NumTicks=5;
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    legendInfo={'t = 0 [s]'; 't = 0.25 [s]';'t = 1 [s]';'Analytical solution at t = 0 [s]'};
    legend(legendInfo,'Interpreter','latex','Location','Best');
    grid
    box on
    %print(strcat(path,'FreeFallingCube_Memory'), '-depsc')
    hold off   



figure(3)
plot(time,mean_Hydrostatic(:,3))
    %axis([0 0 0 0])
    title('')
    set(gca,'FontSize',22);
    xlabel('Time [s]','FontSize',22,'Interpreter','latex');
    ylabel('Hydrostatic Pressure [Pa]','FontSize',22,'Interpreter','latex');
    NumTicks=5;
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    %legendInfo={ 'SPH';'SPH Pressure no init'; 'SPH $\alpha = 0.1$';'SPH $B = 128$';'XSPH $\epsilon = 0.5$'};
    %legend(legendInfo,'Interpreter','latex','Location','Best');
    grid
    box on
    hold off
    
% Save data
DATA.name = dirName(1).name;
DATA.path = nameExperiment;
DATA.nbrFiles = nstep;
DATA.timeStep = timeStep;
DATA.timeWrite = timeWrite;
DATA.timeSimu = timeSimu;
DATA.currTimeSimu = currTimeSimu;
DATA.CPUtime = CPU_Time;
DATA.memory = Memory;
DATA.memoryPeak = Memory_Peak;
DATA.time = time;
DATA.Height = Height;
DATA.H = H;
DATA.meanDensity = mean_Density;
DATA.stdDensity = std_Density;
DATA.meanHydrostatic = mean_Hydrostatic;
DATA.stdHydrostatic = std_Hydrostatic;
DATA.DeltaLog = Delta_log;
DATA.epsilonLog = epsilon_log;
DATA.freq = freq;

save(strcat(nameExperiment,strcat('/',dirName(1).name(1:end-13))), 'DATA')



%% Piston (Gas)
%  ************************************************************************
case 3

    for i=1 : nstep 
        
        % Open File nbr i
        %disp(dirName(i+1).name);
        filename=strcat(nameExperiment,'/',dirName(i).name);
        Experiment = importdata(filename); % Import Data
        
        % Pressure mean of the gas inside the piston
        mean_Pressure(i) = mean(Experiment.data(1:limit(1),8));
        std_Pressure(i)  =  std(Experiment.data(1:limit(1),8));
        index = find(Experiment.data(limit(1):end,6) == -0.05);
        volumePiston(i) = 0.4*0.4*(mean(Experiment.data(index,3)));
    end
transpose(volumePiston)

figure(1)
hold on
    errorbar(time, mean_Pressure, std_Pressure, '*','LineWidth', 2)

    %axis([0 0.2 -0.1e5 4e5])
    set(gca,'FontSize',22);
    xlabel('Time [s]','FontSize',22,'Interpreter','latex');
    ylabel('Mean Pressure [Pa]','FontSize',22,'Interpreter','latex');
    NumTicks=5;
    L = get(gca,'XLim');
    %set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    L = get(gca,'YLim');
    %set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    %legendInfo={'t = 0 [s]'; 't = 1 [s]';'t = 2 [s]';'Analytical solution'};
    %legend(legendInfo,'Interpreter','latex','Location','Best');
    grid
    box on
    %print(strcat(path,'FreeFallingCube_Memory'), '-depsc')
    hold off   

    
    figure(2)
    hold on
    errorbar(volumePiston, mean_Pressure, std_Pressure, '*','LineWidth', 2)
    %axis([0 0.2 -0.1e5 4e5])
    set(gca,'FontSize',22);
    xlabel('Volume $[m^3]$','FontSize',22,'Interpreter','latex');
    ylabel('Mean Pressure [Pa]','FontSize',22,'Interpreter','latex');
    NumTicks=5;
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    %legendInfo={'t = 0 [s]'; 't = 1 [s]';'t = 2 [s]';'Analytical solution'};
    %legend(legendInfo,'Interpreter','latex','Location','Best');
    grid
    box on
    %print(strcat(path,'FreeFallingCube_Memory'), '-depsc')
    hold off 
    
figure(3)
hold on
    plot(time, mean_Pressure.*volumePiston, '*','LineWidth', 2)

    %saxis([0 0.2 -0.1e5 4e5])
    set(gca,'FontSize',22);
    xlabel('Time [s]','FontSize',22,'Interpreter','latex');
    ylabel('P*V $[kg * m^2 / s^2]$','FontSize',22,'Interpreter','latex');
    NumTicks=5;
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    %legendInfo={'t = 0 [s]'; 't = 1 [s]';'t = 2 [s]';'Analytical solution'};
    %legend(legendInfo,'Interpreter','latex','Location','Best');
    grid
    box on
    %print(strcat(path,'FreeFallingCube_Memory'), '-depsc')
    hold off      
% Save data
DATA.name = dirName(1).name;
DATA.path = nameExperiment;
DATA.nbrFiles = nstep;
DATA.timeStep = timeStep;
DATA.timeWrite = timeWrite;
DATA.timeSimu = timeSimu;
DATA.currTimeSimu = currTimeSimu;
DATA.CPUtime = CPU_Time;
DATA.memory = Memory;
DATA.memoryPeak = Memory_Peak;
DATA.time = time;
DATA.meanHydrostatic = mean_Pressure;
DATA.stdHydrostatic = std_Pressure;

save(strcat(nameExperiment,strcat('/',dirName(1).name(1:end-13))), 'DATA')


%% CrashCube
%  ************************************************************************
case 4

    
    % Parameters
    g = 9.81;
    z0_center = 10; %[m]
    
    for i=1 : nstep 
        % Open File nbr i
        %disp(dirName(i+1).name);
        filename=strcat(nameExperiment,'/',dirName(i).name);
        Experiment = importdata(filename); % Import Data
        
        % Compute the Mean in Z coordinate
        MeanZ_experiment(i) = mean(abs(Experiment.data(:,3)));
    end
 
    
    % Plot
    figure(1)
    hold on
    plot(time, MeanZ_experiment);
        title('')
        xlabel('Time [s]')
        ylabel('Z')
        text =  strcat('Simulation (time step = ', num2str(timeStep), ')');
        legend( text)
        grid
        %print('FreeFallingCube_error', '-depsc')
    hold off  
    

    
% Save data
DATA.name = dirName(1).name;
DATA.path = nameExperiment;
DATA.nbrFiles = nstep;
DATA.timeStep = timeStep;
DATA.timeWrite = timeWrite;
DATA.timeSimu = timeSimu;
DATA.currTimeSimu = currTimeSimu;
DATA.CPUtime = CPU_Time;
DATA.memory = Memory;
DATA.memoryPeak = Memory_Peak;
DATA.time = time;
DATA.meanZexperiment = MeanZ_experiment;

save(strcat(nameExperiment,strcat('/',dirName(1).name(1:end-13))), 'DATA')

%% Dam break
%  ************************************************************************
case 5

     L_init = max(InitExperiment.data(1:limit(1),2))- min(InitExperiment.data(1:limit(1),2));
     
     for i=1 : nstep
        
        % Open File nbr i
        %disp(dirName(i+1).name);
        filename=strcat(nameExperiment,'/',dirName(i).name);
        Experiment = importdata(filename); % Import Data
        
        % Compute time of the experiment
        time(i) = (i-1)*timeWrite; % Time 
        
    
    end
    
    
    % Put file in a new directory
    %     if(input('MoveFile to an other folder (1 for yes, 0 for no): ' ))
    %     mkdir ../build/Results/FreeFallingCube
    %     movefile ../build/Results/FreeFallingCube_*.txt ../build/Results/FreeFallingCube
    %     movefile ../build/Results/FreeFallingCube_*.vtk ../build/Results/FreeFallingCube
    %     end

    
    for i=1 : nstep 
        % Open File nbr i
        %disp(dirName(i+1).name);
        filename=strcat(nameExperiment,'/',dirName(i).name);
        Experiment = importdata(filename); % Import Data
        
        L = max(Experiment.data(1:limit(1),2))- min(Experiment.data(1:limit(1),2));
        part = find (Experiment.data(1:limit(1),2) > (min(Experiment.data(1:limit(1),2))+0.9*L) & Experiment.data(1:limit(1),2) <= (min(Experiment.data(1:limit(1),2))+L));  
        v(i) = mean(Experiment.data(part(:),5));
        %part2 = find(Experiment.data(1:limit(1),2)==max(Experiment.data(1:limit(1),2)));
        %v2(i) = mean(Experiment.data(part2(:),5));
        X(i) = L;
        
    end
    
    
    X = X./L_init;
    t = time.*sqrt(2*9.81/L_init);
    
    t1_exp = [0 0.383261 0.773638 1.15393 1.52373 1.93278 2.3185 2.70362 3.08399];
    X1_exp = [1 1.11204 1.25257 1.50677 1.89359 2.24755 2.61080 3.00249 3.62163];
    t2_exp = [0 0.4063 0.8591 1.1862 1.4280 1.6316 1.8201 1.9849 2.1882 2.3149 2.5031 2.6372 2.8178 2.9672 3.0937];
    X2_exp = [1 1.1122 1.2201 1.4359 1.6748 1.8944 2.1044 2.3333 2.5718 2.7814 3.0009  3.2294 3.4394 3.6728 3.8918];
    t3_exp = [0 0.3986 0.8206 1.5622 1.9001 2.2146 2.5602 2.89];
    X3_exp = [1 1.1122 1.2198 1.8939 2.3326 2.7853 3.2241 3.6769];
    
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
    h=legend('Numerical results','Experiment 1','Experiment 2','Experiment 3');
    set(h,'Interpreter','Latex')
    
    figure(2)
    hold on
    plot(time,v)
    grid on
    box on
    ylabel('$u_y$ [m/s]','Interpreter','Latex')
    xlabel('Time [s]','Interpreter','Latex')
    set(gca,'FontSize',25)
    set(h,'Interpreter','Latex')

%% Not Valid Experiment
%  ************************************************************************
otherwise
        disp('Not Valid Experiment')
end

sucess = 0;
end




