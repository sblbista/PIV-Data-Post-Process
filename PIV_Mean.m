clear,clc
format long

% Getting path to datasets
pathToFolder = './AMT00HighReAll';
files = dir( fullfile(pathToFolder,'*.dat') );

% Reading all files
NFiles = numel(files);
dataV = cell(NFiles,1); %set up the length of the dataV

% Initialize the sum matrix
sumU=0.0;
sumV=0.0;

disp('<strong>Reading files...</strong>')
for i=1:numel(files)
    fid = fopen(fullfile(pathToFolder,files(i).name), 'rt');
    H = textscan(fid, '%s', 0, 'Delimiter','\n'); 
    lineformat = repmat('%f',1,4);
    C = textscan(fid, lineformat, 'delimiter', ',','HeaderLines', 1, 'CollectOutput',1);
    fclose(fid);
 
    % Data columns for Uinst  and Vinst
    dataV{i} = [C{1,1}(:,3) C{1,1}(:,4)];
    sumU= sumU + dataV{i}(:,1);
    sumV= sumV + dataV{i}(:,2);
end 

% Getting shape of x and y
shapeX= size(unique(C{1,1}(:,1)));
shapeY= size(unique(C{1,1}(:,2)));

% Getting y data points
y= unique(C{1,1}(:,2));
x= unique(C{1,1}(:,1));

disp('<strong>Calculating PIV statistics...</strong>')

%% Calculating mean velocities for the datasets
Umean=sumU/NFiles;
Vmean=sumV/NFiles;

% Reshaping mean U and V into (216,86) matrix
UMean= reshape(Umean,[shapeX(1),shapeY(1)]);
VMean= reshape(Vmean, [shapeX(1),shapeY(1)]);
disp('<strong>Mean Calculation Completed!</strong>');

%% Calculating velocity fluctuations

% for k= 1:NFiles
%     UPrime= (dataV{k}(:,1)) - Umean(k); 
%     VPrime= (dataV{k}(:,2)) - Vmean(k);
% end
% 
% mean(UPrime)

%% Calculate vorticity
% Initialize dvdx and dudy matrices
dvdx= zeros(size(UMean));
dudy= zeros(size(UMean));

dx= x(2) - x(1);
dy= y(2) - y(1);


for i= 1:numel(x)
    for j= 1:numel(y)
        if i==1
            % Implement Forward Difference Scheme
            dudy(i,j)= (UMean(i+1,j) - UMean(i,j)) / (dy);
            dvdx(i,j)= (VMean(i+1,j) - VMean(i,j)) / (dx);
        elseif i==numel(x)
            % Implement Backward Difference Scheme
            dudy(i,j)= (UMean(i,j) - UMean(i-1,j)) / (dy);
            dvdx(i,j)= (VMean(i,j) - VMean(i-1,j)) / (dx);
        else
            % Central differnce method
            dudy(i,j)= (UMean(i+1,j) - UMean(i-1,j)) / (2 * dy);
            dvdx(i,j)= (VMean(i+1,j) - VMean(i-1,j)) / (2 * dx);
        end
    end
end

% Vorticity
w= dvdx - dudy;

w= reshape(w,size(Umean));

%% Export to TECPLOT format 
disp('Working on results files...');
PIVStats = [C{1,1}(:,1) C{1,1}(:,2) Umean Vmean w];
filename = 'PIV_Results.dat';
fid = fopen(filename, 'w');
fprintf(fid, 'TITLE=%s\n', filename);
fprintf(fid, 'VARIABLES= X, Y, U, V Wz \n');
fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', shapeX(1), shapeY(1));
dlmwrite(filename, PIVStats, '-append', 'delimiter', ' ');
fclose(fid);

disp('<strong>PIV statistics completed</strong>');





