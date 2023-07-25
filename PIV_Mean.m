%% This code calculates the turbulent quantities from experimental data %% 
clear,clc
format long

% Getting path to datasets
pathToFolder = './Slot_Output'; % Enter path to folder
files = dir( fullfile(pathToFolder,'*.dat') );

% Reading all files
NFiles = numel(files);
dataV = cell(NFiles,1); %set up the length of the dataV

% Initialize the sum matrix
sumU= 0.0;
sumV= 0.0;

disp('<strong>Reading files...</strong>')
%% Read data and calculate the sum
for i=1:numel(files)
    fid = fopen(fullfile(pathToFolder,files(i).name), 'rt');
    H = textscan(fid, '%s', 5, 'Delimiter','\n');  % Header lines to skip
    lineformat = repmat('%f',1,5); % Number of unique columns
    C = textscan(fid, lineformat, 'delimiter', ',','HeaderLines', 1, 'CollectOutput',1);
    fclose(fid);
 
    % Data columns for Uinst  and Vinst
    dataV{i} = [C{1,1}(:,3) C{1,1}(:,4)];
    sumU= sumU + dataV{i}(:,1);
    sumV= sumV + dataV{i}(:,2);
    u= dataV{i}(:,1);
    v= dataV{i}(:,2);
end 
% Getting y data points
y= unique(C{1,1}(:,2));
x= unique(C{1,1}(:,1));

% Getting shape of x and y
rows= size(x,1);
cols= size(y,1);

disp('<strong>Calculating PIV statistics...</strong>')

%% Calculating mean velocities for the datasets
Umean= sumU/NFiles;
Vmean= sumV/NFiles;

% Reshaping mean U and V into (x,y) matrix
UMean= reshape(Umean,[rows,cols]);
VMean= reshape(Vmean, [rows,cols]);
disp('<strong>MEAN VELOCITY CALCULATION COMPLETE!</strong>');
x_y_size= size(UMean);
tecSize= size(Umean);

%% Velocity fluctuations for U and V
U_flucs_sum = 0;%[]; % This needs to be an array of fluctuations
V_flucs_sum = 0;%[]; 

% Initializing the Urms and Vrms
Urms = 0;
Vrms = 0;

% Initialzing the UVres
UVres = 0;

for k=1:NFiles
   % Calculating the velocity fluctuations for individual images
    U_flucs_sum= (dataV{k}(:,1)-Umean) + U_flucs_sum;
    V_flucs_sum= (dataV{k}(:,2)-Vmean) + V_flucs_sum;

   % calculating the root mean squares of U, V
   Urms =(dataV{k}(:,1)-Umean).^2 + Urms;
   Vrms =(dataV{k}(:,2)-Vmean).^2 + Vrms;

   % calculating the reynold stress, UV
   UVres= (dataV{k}(:,1)-Umean).*(dataV{k}(:,2)-Vmean)+ UVres;
end    

% Calculating Urms and Vrms
Urms = sqrt(Urms/NFiles);
Vrms = sqrt(Vrms/NFiles);
disp('<strong>RMS CALCULATION COMPLETE!</strong>');

% Calculating UVres
UVres = UVres/NFiles;


% calculating the velocity fluctuations for all images
Uprime = U_flucs_sum/NFiles;
Vprime = V_flucs_sum/NFiles;

% Normal Reynolds Stress
uuStress= (Urms)/NFiles;
vvStress= (Vrms)/NFiles;
uvStress= UVres/NFiles;
disp('<strong>REYNOLDS STRESS CALCULATION COMPLETE!</strong>');

% calculate Turbulent Kinetic Energy
TKE= 0.5 * (uuStress + vvStress);



%% Calculate vorticity
% Initialize dvdx and dudy matrices
u=reshape(u,[rows,cols]);
v=reshape(v,[rows,cols]);

dx= x(2) - x(1);
dy= y(2) - y(1);

%% Calculate the velocity gradient tensor and Swirling Strength 
dvdx= zeros(size(x_y_size));
dudy= zeros(size(x_y_size));
dvdy= zeros(size(x_y_size));
dudx= zeros(size(x_y_size));

[dudx, dudy]= gradient(UMean,dx,dy); % Change UMean to u (instantenous velocity)
[dvdx, dvdy]= gradient(VMean,dx,dy); % Change VMean to v (instantenous velocity)


% Calculate Turbulent Production
Prod= - (reshape(uuStress,[rows,cols]).*dudx + reshape(uvStress,[rows,cols]).*dvdx + reshape(uvStress,[rows,cols]).*dudy + reshape(vvStress,[rows,cols]).*dvdy);
Prod= reshape(Prod,tecSize);
% Calculate Vorticity
wz= (dvdx - dudy);

for i= 1:numel(x)
    for j= 1:numel(y)
        % Velocity Gradient Tesor
        D2D= [dudx(i,j) dudy(i,j); dvdx(i,j) dvdy(i,j)];
        temp_eig= eig(D2D);
        lambda(i,j)= imag(temp_eig(1,1));
        % Swirling Strength
        lambda1(i,j)=lambda(i,j)*sign(wz(i,j));
        P(i,j)= trace(D2D); % 0 for incompressible flow 
        % Q-Criterion
        Q(i,j)= -(dudy(i,j).*dvdx(i,j)) + (dudx(i,j).*dvdy(i,j));
    end
end

disp('<strong>VORTICITY CALCULATION COMPLETE!</strong>')

lambda1= lambda1';
Wz= reshape(wz,tecSize);
lambda1R= reshape(lambda1,tecSize);

tec= input('Write data to TECPLOT? [0/1]');
if tec==1
    %% Export to TECPLOT format 
    disp('Working on results files...');
    PIVStats = [C{1,1}(:,1) C{1,1}(:,2) Umean Vmean uuStress vvStress uvStress TKE Prod Wz lambda1R];
    filename = 'Slot0.dat';
    fid = fopen(filename, 'w');
    fprintf(fid, 'TITLE=%s\n', filename);
    fprintf(fid, "VARIABLES= X, Y, U, V U'U' V'V' U'V' TKE PRODUCTION VORT LAMBDA_CI \n");
    fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', rows, cols);
    dlmwrite(filename, PIVStats, '-append', 'delimiter', ' ');
    fclose(fid);
    
    disp('<strong>PIV statistics completed!</strong>');
elseif tec==0
end

% Write the variables to a file
flag= input('Write data as CSV file? [0/1]');

if flag==1
    %% Write all the variables as CSV file
    M= [C{1,1}(:,1) C{1,1}(:,2) Umean Vmean Urms Vrms UVres];
    csvwrite('outputPIV.csv',M);
end

% Perform double averging of the mean quantities
if flag==1
    [sm1 sm2]=size(M(:,1));
    [sm3 sm4]=size(M(1,:));
    for i=1:sm4
        M1=reshape(M(:,i),[r2,r1]);
        for j=1:r2
            M2(j,i)=mean(M1(j,:));
        end
    end
end


    

