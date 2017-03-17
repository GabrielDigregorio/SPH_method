%**************************************************************************
% Plot Playground:
%**************************************************************************
function sucess = Playground(filename, n, path)


close all;
set(groot,'defaultLineLineWidth',2)

% Check path (optional arument)
if nargin < 2
   path = '';
end

% Import data at t=0
i=0;
filename=strcat('../build/Results/', filename,'_',num2str(0,'%08i'),'.txt')
InitExperiment = importdata(filename);

        Str1 = char(InitExperiment.textdata(11));
             Key1 = 'Domain (lower l) : ';
        Str2 = char(InitExperiment.textdata(12));
             Key2 = 'Domain (upper u) : ';
        Str3 = char(InitExperiment.textdata(13));
             Key3 = 'Number of Particules (nFree/nMoving/nFixed) : ';             
        Index1 = strfind(Str1, Key1); Index2 = strfind(Str2, Key2); Index3 = strfind(Str3, Key3);
        l     = sscanf(Str1(Index1(1) + length(Key1):end), '%g %g %g', 3);
        u     = sscanf(Str2(Index2(1) + length(Key2):end), '%g %g %g', 3);
        limit = sscanf(Str3(Index3(1) + length(Key3):end), '%g %g %g', 3);

%% SEE GEOMETRY

    disp(['What would you like to plot: ']); % Display Possibilities
    disp(['1) All particules (Free, Moving and Fixed)']);
    disp(['2) Free particules']);
    disp(['3) Free particules and Moving particules']);disp([' ']);
    n = input('Choice: ');
        
switch n

    case 1
        scatter3(InitExperiment.data(:,1),InitExperiment.data(:,2),InitExperiment.data(:,3),10,InitExperiment.data(:,3))
    case 2
        scatter3(InitExperiment.data(1:limit(1),1),InitExperiment.data(1:limit(1),2),InitExperiment.data(1:limit(1),3),10,InitExperiment.data(1:limit(1),3))
    case 3
        scatter3(InitExperiment.data(1:limit(1)+limit(2),1),InitExperiment.data(1:limit(1)+limit(2),2),InitExperiment.data(1:limit(1)+limit(2),3),10,InitExperiment.data(1:limit(1)+limit(2),3))
end

axis([l(1) u(1) l(2) u(2) l(3) u(3)])
colormap(jet);
colorbar;



sucess = 0;
end






