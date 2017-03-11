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

% Choices
disp(['Would you like to:']);
disp(['    1) Compile the project; ']);
disp(['    2) Run an experiment; ']);
disp(['    3) Analyse an experiment; ', ' ']);
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
        if(strcmp(name,'DESKTOP-31TT348'))
           p = input('Parameter file name: ','s'); p = strcat(p, '.txt ')
           g = input('Geometry file name: ','s'); g = strcat(g, '.txt ')
           system(char(strcat({'"../build/sph.exe"'},{' '},p,{' '},g,{' '}, {'FreeFallingCube'})))
        else disp(['You are not allowed to launch an experiment... ']);
        end

    case 3 % Analyse:
        disp(['Experiments: ']); % Display Possibilities
        disp(['1) Free Falling Cube']);
        disp(['2) Free Falling Cube with random particules']);
        disp(['3) Stationary Tank']);disp([' ']);
        n = input('Enter the number of the experiment: ');
        exit = TestCases(n, path)
        
    otherwise disp(['Not valid choice']);
end






