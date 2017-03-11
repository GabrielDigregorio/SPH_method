%**************************************************************************
% Test Cases:
%                 - Free falling cube under gravity
%                 - Stationary tank
%                 - Crashed Cube
%**************************************************************************
clc; clear all; close all;
set(groot,'defaultLineLineWidth',2)

% Get Environement:
name = getenv('COMPUTERNAME');
if (strcmp(name,'DESKTOP-31TT348')) 
     path = 'C:\Users\gabri\Dropbox\Applications\ShareLaTeX\SPH_PROJECT\';
     disp(['Welcome GabyGab: ']);
else path = '';
end

% Display Possibilities
disp(['Experiments: ']);
disp(['1) Free Falling Cube']);
disp(['2) Free Falling Cube with random particules']);
disp(['3) Stationary Tank']);disp([' ']);
n = input('Enter the number of the experiment: ');

switch n


%% Free Falling cube
%  ************************************************************************
case 1

    % Parameters
          g = 9.81;
%         z0_center = 100; %[m]
          nstep = 1; %[-]
          timeStep = 0.1; %[s]
    % Cube:
%         L=10;
%         W=10;
%         H=10;
%         r=0;

    % Import data at t=0
    filename=strcat('../Results/FreeFallingCube_',num2str(0,'%08i'),'.txt')
    InitExperiment = importdata(filename);
        
    for i=1 : nstep 
        
        % Open File nbr i
        filename=strcat('../Results/FreeFallingCube_',num2str(i,'%08i'),'.txt')
        Experiment = importdata(filename); % Import Data
        
        % Compute the error to analytical solution
        time(i) = (i-1)*timeStep; % Time 
        Analytic = InitExperiment.data(:,3) - g*(time(i)^2)/2; % MRUA
        error(i) = mean(abs(((Experiment.data(:,3)-Analytic)./Analytic)))*100; % error [%]
        XY_move(i) = sqrt(mean(  (InitExperiment.data(:,1) - Experiment.data(:,1)).^2  ...
                              +  (InitExperiment.data(:,2) - Experiment.data(:,2)).^2   ));
    
        % Get memory consumption
        Str1 = char(Experiment.textdata(7)); Key1 = 'Memory Usage :';
        Str2 = char(Experiment.textdata(8)); Key2 = 'Peak :';
        Index1 = strfind(Str1, Key1); Index2 = strfind(Str2, Key2);
        Memory(i)      = sscanf(Str1(Index1(1) + length(Key1):end), '%g', 1)
        Memory_Peak(i) = sscanf(Str2(Index2(1) + length(Key2):end), '%g', 1)
    
    end

    % Plot
    figure(1)
    hold on
    plot(time, error, '*');
    % polyfit...
        %axis([0 0 0 0])
        title('')
        xlabel('Time [s]')
        ylabel('Error in %')
        legend('')
        grid
        print('FreeFallingCube_error', '-depsc')
    hold off  
    
    figure(2)
    hold on
    plot(time, XY_move, '*');
        %axis([0 0 0 0])
        title('')
        xlabel('Time [s]')
        ylabel('Error in %')
        legend('')
        grid
        %print(strcat(path,'FreeFallingCube_XY_move'), '-depsc')
    hold off  




%% Free Falling cube Random
%  ************************************************************************
case 2

    % Parameters
    g = 9.81;
    z0_center = 100; %[m]
    nstep = 1; %[-]
    timeStep = 0.1; %[s]
    % Cube:
    L=10;
    W=10;
    H=10;
    r=20;


    % Import data at t=0
    filename=strcat('../Results/FreeFallingCube_',num2str(0,'%010i'),'.txt')
    InitExperiment = importdata(filename);
        
    for i=1 : nstep 
        filename=strcat('../Results/FreeFallingCube_',num2str(i,'%010i'),'.txt')
        Experiment(i) = importdata(filename);
        time(i) = (i-1)*timeStep;
        Analytic = InitExperiment.data(:,3) - g*(time(i)^2)/2; % MRUA
        error(i) = mean(abs(((Experiment(i).data(:,3)-Analytic)./Analytic)))*100; % [%]
        XY_move(i) = sqrt(mean(  (InitExperiment.data(:,1) - Experiment(i).data(:,1)).^2  ...
                              +  (InitExperiment.data(:,2) - Experiment(i).data(:,2)).^2   ));
    end

    % Plot
    figure(10)
    hold on
    plot(time, error, '*');
        %axis([0 0 0 0])
        title('')
        xlabel('Time [s]')
        ylabel('Error in %')
        legend('')
        grid
        %print(strcat(path,'FreeFallingCubeRandom'), '-depsc')
    hold off  

    
    

%% Stationary tank
%  ************************************************************************
case 3

    % Parameters

    % Cube

    % Import data


    % Plot
    figure(20)
    hold on
    plot();
        %axis([0 0 0 0])
        title('')
        xlabel('')
        ylabel('')
        legend('')
        grid
        %print(strcat(path,'StationaryTank'), '-depsc')
    hold off  

    
    
%% Not Valid Experiment
%  ************************************************************************
otherwise
        disp('Not Valid Experiment')
end




