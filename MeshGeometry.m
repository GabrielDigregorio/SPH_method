%clear all; clc; close all;
set(groot,'defaultLineLineWidth',2)

%% SEE GEOMETRY
A = importdata('Playground.txt');

scatter3(A(:,1),A(:,2),A(:,3),10,A(:,3))
axis([-20 20 -20 20 -20 20])
colormap(jet);
colorbar;


%% ANALYS TIME ALGO SEARCH
% A = importdata('neighborAnalysis.txt');
% 
% plot(A(:,1),A(:,2),A(:,1),A(:,3))
%         %axis([-216 0 0 2e-3])
%         title('Time Search Algo')
%         xlabel('number of particles [-]')
%         ylabel('time [s]')
%         legend('AllPair','linked-list')
%         grid








