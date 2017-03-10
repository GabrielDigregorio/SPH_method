%**************************************************************************
% Test Cases:
%                 - Free falling cube under gravity
%                 - Stationary tank
%                 - Crashed Cube
%**************************************************************************
clc; clear all; close all;
set(groot,'defaultLineLineWidth',2)

%% Free Falling cube
%  ************************************************************************

% Parameters
g = 9.81;
z0_center = 100; %[m]
nstep = 1; %[-]
timeStep = 0.1; %[s]
% Cube:
L=10;
W=10;
H=10;
r=0;


% Import data
filename=strcat('../Results/FreeFallingCube_0.txt')
InitExperiment = importdata(filename);
    
for i=1 : nstep 
    disp('ok\n')
    filename=strcat('../Results/FreeFallingCube_',num2str(i),'.txt')
    Experiment(i) = importdata(filename);
    time(i) = (i-1)*timeStep;
    Analytic = InitExperiment.data(:,3) - g*(time(i)^2)/2;
    error(i) = mean(abs(((Experiment(i).data(:,3)-Analytic)./Analytic)))*100; % [%]
end

% Plot
figure(1)
hold on
plot(time, error, '*');
    %axis([0 0 0 0])
    title('')
    xlabel('Time [s]')
    ylabel('Error in %')
    legend('')
    grid
    print('FreeFallingCube', '-depsc')
hold off  


%% Stationary tank
%  ************************************************************************
% Parameters

% Cube

% Import data


% Plot
figure(2)
hold on
plot();
    %axis([0 0 0 0])
    title('')
    xlabel('')
    ylabel('')
    legend('')
    grid
    print('StationaryTank', '-depsc')
hold off  







