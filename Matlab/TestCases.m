%**************************************************************************
% Test Cases:
%                 - Free falling cube under gravity
%                 - Stationary tank
%                 - Crashed Cube
%**************************************************************************
function sucess = TestCases(nameExperiment, n, path)

close all;
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

    % Domain and nbr particules
    Str1 = char(InitExperiment.textdata(12));
    Key1 = 'Domain (lower l) : ';
    Str2 = char(InitExperiment.textdata(13));
    Key2 = 'Domain (upper u) : ';
    Str3 = char(InitExperiment.textdata(14));
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
    
    %     figure(3)
    %     hold on
    %     plot(time, XY_move, '*');
    %         %axis([0 0 0 0])
    %         title('')
    %         xlabel('Time [s]')
    %         ylabel('Error in %')
    %         legend('')
    %         grid
    %         %print(strcat(path,'FreeFallingCube_XY_move'), '-depsc')
    %     hold off  
    
    figure(4)
    hold on
    plot(time, Memory*1.25e-7, time, Memory_Peak*1.25e-7);
        title('Memory Consumption')
        xlabel('Time [s]')
        ylabel('Memory Consumption [MB]')
        legend('RSS Memory', 'RSS Memory Peak')
        grid
        %print(strcat(path,'FreeFallingCube_Memory'), '-depsc')
    hold off 
    
    figure(5)
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
    nbrWindows = 4;
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
        for j=1:nbrWindows+1
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
    for i=[floor(length(H(:,1))/2)  length(H(:,1))]
        errorbar(H(i,:), mean_Hydrostatic(i,:), std_Hydrostatic(i,:), '*','LineWidth', 2)
    end
    plot(H(i,:), 1000*9.81*(height_max -H(i,:)), 'color', 'k')
    %axis([0 1 2 3])
    set(gca,'FontSize',22);
    xlabel('Height [m]','FontSize',22,'Interpreter','latex');
    ylabel('Hydrostatic Pressure [Pa]','FontSize',22,'Interpreter','latex');
    NumTicks=5;
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    legendInfo={'t = 0 [s]'; 't = 1 [s]';'t = 2 [s]';'Analytical solution'};
    legend(legendInfo,'Interpreter','latex','Location','Best');
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
        volumePiston(i) = 0.38*0.38*(mean(Experiment.data(index,3))-0.1);
    end


figure(1)
hold on
    errorbar(time, mean_Pressure, std_Pressure, '*','LineWidth', 2)

    axis([0 0.2 -0.1e5 4e5])
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
DATA.CPUtime = CPU_Time;
DATA.memory = Memory;
DATA.memoryPeak = Memory_Peak;
DATA.time = time;
DATA.meanZexperiment = MeanZ_experiment;

save(strcat(nameExperiment,strcat('/',dirName(1).name(1:end-13))), 'DATA')

%% Not Valid Experiment
%  ************************************************************************
otherwise
        disp('Not Valid Experiment')
end

sucess = 0;
end




