%**************************************************************************
% GUI: Pseudo Graphical Interface
%                 - Compile the project
%                 - Run an experiment
%                 - Analyse an experiment
%**************************************************************************
clc; clear all; close all;

% Get Environement:
name = getenv('COMPUTERNAME');
if (strcmp(name,'DESKTOP-31TT348')) 
     path = 'C:\Users\gabri\Dropbox\Applications\ShareLaTeX\SPH_PROJECT\';
     disp(['Welcome GabyGab: ']); disp(['*****************']);disp([' ']);
else path = ''; 
end

loop=1;

while (loop==1)
% Choices
disp(['Would you like to:']);
disp(['    1) Compile the project; ']);
disp(['    2) Run an experiment; ']);
disp(['    3) Plot the Playground of an experiment; ']);
disp(['    4) Analyse an experiment; ']);
disp(['    5) Edit a new experiment; ']);
disp(['    6) Exit; ', ' ']);
choice = input('Enter your choice number: ' ); disp(' ');

switch(choice)
    case 1 % Compile:
        cd ..; cd build\;
        if(ispc) system('mingw32-make')
        elseif(isunix) system('make')
        else disp(['You are not allowed to compile... ']);
        end
        cd ..; cd Matlab\;
        
    case 2 % Run:
        %if(strcmp(name,'DESKTOP-31TT348'))
           disp(['Choose among the list of Playgrounds: ']);
           cd ..\Playgrounds\
           
           disp(['Clic on Parameter file: ']);
           [p, pathname] = uigetfile({'*.kzr'},'\Playgrounds\');
           
           disp(['Clic on Geometry file: ']);
           [g, pathname] = uigetfile({'*.kzr'},'\Playgrounds\');
           
           p = strcat('..\Playgrounds\',p)
           g = strcat('..\Playgrounds\',g)
           
           cd ../build/Results
           nameExperiment = input('Enter the name of the experiment: ','s');
           mkdir(nameExperiment);
           nameExperiment = strcat(nameExperiment,'/', nameExperiment);
           cd ..
           system(char(strcat({'"sph.exe"'},{' '},p,{' '},g,{' '}, {nameExperiment})))
           cd ..\Matlab\;
        %else disp(['You are not allowed to launch an experiment... ']);
        %end
        
    case 3 % Playground:
        disp(['Choose an experiment : ']);
        nameExperiment = uigetdir('../build/Results/');
        Playground(nameExperiment, 1, path);
        
    case 4 % Analyse:
        disp(['Experiments Types: ']); % Display Possibilities
        disp(['1) Free Falling Cube']);
        disp(['2) Swimming Pool']);
        disp(['3) not defined']);disp([' ']);
        n = input('Enter the TYPE of the experiment: ');disp([' ']);
        disp(['Clic on directory: ']);
        nameExperiment = uigetdir('../build/Results/');
        exit = TestCases(nameExperiment, n, path);
    
    case 5 % Edit File
        cd ..; cd Playgrounds\;
        p = input('NEW Parameter file name: ','s'); p = strcat(p, '.kzr ');
        if (exist(p, 'file')==0) 
            fid = fopen( p, 'wt' );
            edit(p);
        end
        g = input('NEW Geometry file name: ','s'); g = strcat(g, '.kzr ');
        if (exist(g, 'file')==0) 
            fid = fopen( g, 'wt' );
            edit(g);
        end
        input('Enter: '); fclose(p); fclose(g);
        cd ..; cd Matlab\;
        
    case 6 % Exit
        loop=0;
    otherwise disp(['Not valid choice']); disp(' ');
end

end






