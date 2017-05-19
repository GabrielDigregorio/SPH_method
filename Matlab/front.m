function [time,t,X,t1_exp,X1_exp,t2_exp,X2_exp,t3_exp,X3_exp,v] = front(nameExperiment)


set(groot,'defaultLineLineWidth',2)

%addpath(genpath('../build/Results/'))

    % Number of files to read
     dirName = dir([nameExperiment, '/*.txt']); %list all directory with.txt
    nstep = length(dir([nameExperiment, '/*.txt']))-1; %[-]

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
    L_init = max(InitExperiment.data(1:limit(1),2))- min(InitExperiment.data(1:limit(1),2));
    % figure(100)
    %     plot(InitExperiment.data(:,3), InitExperiment.data(:,8))

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
    
end

