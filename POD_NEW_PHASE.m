%Developed by Sabal Bista (bistas@uwindsor.ca)
% ---This code reads PIV data and computes the POD modes, energy distrubution,
% and mode coefficients using SVD technique. It also reconstructs the
% velocity field using specified energy criteria and displays the velocity
% magnitude. It also computes the phase angles based on the first two
% dominant modes.---%  
% Feel free to modify the code to your liking % 

clc; clear; close all;

pathToFolder = 'HCG100_ALL/HCG100ZA_COMBINEDFOV';
files = dir( fullfile(pathToFolder,'*.dat') );
NFiles =  500; %numel(files);
%read all files
dataV = cell(NFiles,1); %set up the length of the dataV

sumU= 0.0;
sumV= 0.0;
% Initialize the sum matrix
%% READ THE AND ARRANGE THEM IN X x Yx T 3D ARRAY 
disp('<strong>Reading files...</strong>')
for i=1:NFiles
    fid = fopen(fullfile(pathToFolder,files(i).name), 'rt');
    H = textscan(fid, '%s', 2, 'Delimiter','\n');  % Header lines to skip
    lineformat = repmat('%f',1,4); % Number of unique columns
    C = textscan(fid, lineformat, 'delimiter', ',','HeaderLines', 1, 'CollectOutput',1);
    fclose(fid);
    
    % Data columns for Uinst  and Vinst
    dataV{i} = C{1,1};
    x= dataV{i}(:,1);
    y= dataV{i}(:,2);
    dataV{i} = [C{1,1}(:,3) C{1,1}(:,4)];
    sumU= sumU + dataV{i}(:,1);
    sumV= sumV + dataV{i}(:,2);
    disp(['Time instance read: ' num2str(i) ]);
end

Umean= sumU/NFiles;
Vmean= sumV/NFiles;
disp('<strong>MEAN COMPUTED!</strong>')


yN= unique(C{1,1}(:,2));
xN= unique(C{1,1}(:,1));
rows= size(xN,1);
cols= size(yN,1);

% Remove the mean from the data
k= 0;
for j= 1:NFiles
    k= k+1;
    Vel_i= dataV{k}(:,1);
    Vel_j= dataV{k}(:,2);
    ufluc= dataV{k}(:,1) - Umean;
    vfluc= dataV{k}(:,2) - Vmean;
    
    UFLUC= reshape(ufluc,[rows,cols]);
    VFLUC= reshape(vfluc,[rows,cols]);
    u_inst=reshape(Vel_i,[rows,cols]);
    v_inst=reshape(Vel_j,[rows,cols]);

    U(:,:,k)= cat(3,UFLUC); %  U(:,i)= cat(3,ufluc) gives (rowsxcols,nFiles) % Don't need to do U_R reshape
    V(:,:,k)= cat(3,VFLUC); %  V(:,i)= cat(3,vfluc)

        
    u(:,:,k)= cat(3,u_inst);
    v(:,:,k)= cat(3,v_inst);
    disp(['Calculating the fluctuating component of time instance: ' num2str(j) ]);
end




disp('<strong>Reshaping Data...</strong>')

X= reshape(x,[rows,cols]);
Y= reshape(y,[rows,cols]);

% a= 0.01; b= 0.05;
% x_indices = find(xN >= a & xN <= b);
% y_indices= find(yN>=0 & yN<=0.00635);

grid_points= rows*cols;
U_R= reshape(U,rows*cols,size(U,3));
V_R= reshape(V,rows*cols,size(V,3));

%% Calculate Reynolds Stress
uu= U_R.^2; uu= mean(uu,2);
vv= V_R.^2; vv= mean(vv,2);
uv= U_R.*V_R; uv= mean(uv,2);

vel_data= (1/(NFiles-1))*[U_R;V_R]; 

% uGap= U(x_indices,y_indices,:);
% vGap= V(x_indices,y_indices,:);
% 
% uGap= reshape(uGap,size(x_indices,1)*size(y_indices,1),size(U,3));
% vGap= reshape(vGap,size(x_indices,1)*size(y_indices,1),size(U,3));

% vel_data= [uGap;vGap];

U_INST= reshape(u,rows*cols,size(U,3));
V_INST= reshape(v,rows*cols,size(V,3));
velocity_data= [U_INST;V_INST];

disp('<strong>COMPUTING MODAL COEFFICIENTS...</strong>')

%% POD ALGORTIHM
[psi,sigma,phi]= svd(vel_data,'econ');
mode_energy= diag(sigma).^2;


%% ENERGY THRESHOLD
energy_threshold = 0.50; % Adjust as needed
cumulative_energy = cumsum(mode_energy) / sum(mode_energy);
num_modes_to_retain = find(cumulative_energy >= energy_threshold, 1);

POD_modes = psi; % (:, 1:num_modes_to_retain)

% M1= reshape(POD_modes(1:1800,1),[75,24]); M2= reshape(POD_modes(1:1800,2),[75,24]); 
% M3= reshape(POD_modes(1:1800,3),[75,24]);M4= reshape(POD_modes(1:1800,4),[75,24]);


PODU= POD_modes(1:grid_points,:);
PODV= POD_modes(grid_points+1:(2*grid_points),:);
POD_MAG= sqrt(PODU.^2 + PODV.^2);

figure;
plot(cumulative_energy(1:10)*100, 'ko--');
xlabel('POD Modes');
ylabel('Cumulative Energy(%)');
% ylim([10,100])

percent_energy(1:10)= (mode_energy(1:10)/(sum(mode_energy)))*100;
figure(2);
plot(percent_energy(1:10),'square -')
xlabel('POD Modes'); ylabel('POD Energy(%)')
title('Energy distributions among the first 10 POD modes')
% Visualize the first few modes 
visMod= input('Visualize modes? [0/1]');
if visMod== 1
    nModes= input('Enter number of modes to visualize(>3): ');
    p= ceil(nModes/2);
    q= ceil(nModes/3);
        
        for i = 1:nModes
            figure(2);
            % hold(q,p,i),'on');
            subplot(q,p,i); 
            contourf((reshape(PODU(1:[rows*cols], i), [rows, cols])'), 20, 'LineColor', 'none');
            title(['POD Mode ', num2str(i)]);
            colorbar; colormap('jet')
        end
end
   
A= vel_data'*POD_modes;
a1= A(:,1); a2= A(:,2);

tec= input('Write data to TECPLOT? [0/1]');
if tec==1
    %% Export to TECPLOT format 
    disp('Working on results files...');
    PODStats = [x y Umean Vmean PODU(:,1) PODU(:,2) PODU(:,3) PODU(:,4) PODU(:,5) PODV(:,1) PODV(:,2) PODV(:,3) PODV(:,4) PODV(:,5)];
    fName= input('Enter file name to export: ',"s");
    filename = [fName,'.dat'];
    fid = fopen(filename, 'w');
    fprintf(fid, 'TITLE=%s\n', filename);
    fprintf(fid, "VARIABLES= X, Y, U, V MODE1U MODE2U MODE3U MODE4U MODE5U MODE1V MODE2V MODE3V MODE4V MODE5V \n");
    fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', rows, cols);
    dlmwrite(filename, PODStats, '-append', 'delimiter', ' ');
    fclose(fid);
    
    disp('<strong>POD MODES EXPORTED!</strong>');
elseif tec==0
end

%% RECONSTRUCT THE DATA
r= 2; %num_modes_to_retain; % reconstruct using r mode
U_red= PODU(:,1:r) * sigma(1:r,1:r) * phi(:,1:r)'; %psi(:,1:r) * sigma(1:r,1:r) * phi(:,1:r)';
V_red= PODV(:,1:r) * sigma(1:r,1:r) * phi(:,1:r)'; 
t= 5; % no. instant
U_REDUCED= reshape(Umean,[rows,cols]) + reshape(U_red(:,t),[rows,cols]);
V_REDUCED= reshape(Vmean,[rows,cols]) + reshape(V_red(:,t),[rows,cols]);
U_REDUCED_MAG= sqrt(U_REDUCED.^2 + V_REDUCED.^2);
contourf(X,Y,U_REDUCED_MAG, 20, 'LineColor','none');
title(['Reconstructed Field using ',num2str(r),' modes'])
colormap('jet')
c= colorbar; c.Label.String= 'U Mag';
axis equal

% EXPORT RECONSTRUCTED DATA
tec= input('Write reconstructed data to TECPLOT? [0/1]');
if tec==1
    RECSTATS= [x y reshape(U_REDUCED,[rows*cols,1]) reshape(V_REDUCED,[rows*cols,1])];
    fName2= input('Enter file name to export: ',"s");
    filename= [fName2,'.dat'];
    fid = fopen(filename, 'w');
    fprintf(fid, 'TITLE=%s\n', filename);
    fprintf(fid, "VARIABLES= X, Y, U_RED, V_RED \n");
    fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', rows, cols);
    dlmwrite(filename, RECSTATS, '-append', 'delimiter', ' ');
    
    disp('<strong>Reconstructed velocity field exported</strong>');
    elseif tec==0
end


%% Phase averaging
% Use the temporal coefficients of the first two modes for phase identification
% Select number of modes to use for phase sorting
n_modes = 2;  % Adjust as needed
lambda1= mode_energy(1)/(sum(mode_energy)); lambda2= mode_energy(2)/(sum(mode_energy));
phase = atan2(sqrt(lambda1)*a2, sqrt(lambda2)*a1);


h1= plot(a1(1:end)/sqrt(2*lambda1), a2(1:end)/sqrt(2*lambda2), 'ko', 'markersize',5);
set(h1, 'markerfacecolor', 'r');
axis equal
xlabel('$a_1(t)/\sqrt{2 \lambda_1}$','Interpreter','latex','FontSize',16)
ylabel('$a_2(t)/\sqrt{2 \lambda_2}$','Interpreter','latex','FontSize',16)

% Normalize phase to [0, 2π]
%phase = mod(phase, 2*pi);
phasedeg= rad2deg(phase);

% Define number of phase bins
n_bins = 16;

% Create phase bins
bin_edges = linspace(0, 2*pi, n_bins+1);

% Sort velocity fields into phase bins
[~, bin_indices] = histc(phase, bin_edges);

% Initialize cell array to store sorted velocity fields
sorted_fields = cell(1, n_bins);



% Sort velocity fields into bins
for i = 1:n_bins
    bin_mask = (bin_indices == i);
    sorted_fields{i} = velocity_data(:, bin_mask);
end

% Compute phase-averaged velocity fields
phase_averaged_fields = cellfun(@(x) mean(x, 2), sorted_fields, 'UniformOutput', false);

% % Reshape phase-averaged fields back to 2D spatial format
% phase_averaged_u = cellfun(@(x) reshape(x(1:rows*cols), [rows, cols]), phase_averaged_fields, 'UniformOutput', false);
% phase_averaged_v = cellfun(@(x) reshape(x(rows*cols+1:end), [rows,cols]), phase_averaged_fields, 'UniformOutput', false);

% Reshape phase-averaged fields back to 2D spatial format
phase_averaged_u = cellfun(@(x) reshape(x(1:rows*cols), [rows*cols,1]), phase_averaged_fields, 'UniformOutput', false);
phase_averaged_v = cellfun(@(x) reshape(x(rows*cols+1:end), [rows*cols,1]), phase_averaged_fields, 'UniformOutput', false);


%% Phase averaged Reynolds stress
u_phase= cell2mat(phase_averaged_u);
v_phase= cell2mat(phase_averaged_v);
u_phase_mean= mean(u_phase,2);
v_phase_mean= mean(v_phase,2);
uu_phase= (u_phase - u_phase_mean).^2; uu_phase= mean(uu_phase,2);
vv_phase= (v_phase - v_phase_mean).^2; vv_phase= mean(vv_phase,2);
uv_phase= ((u_phase - u_phase_mean).*(v_phase - v_phase_mean)); uv_phase= mean(uv_phase,2);

% EXPORT PHASE-AVERAGED DATA
tec= input('Write Reynolds stress data to TECPLOT? [0/1]');
if tec==1
    RECSTATS= [x y uu uv uu_phase uv_phase];
    fName3= input('Enter file name to export: ',"s");
    filename= [fName3,'.dat'];
    fid = fopen(filename, 'w');
    fprintf(fid, 'TITLE=%s\n', filename);
    fprintf(fid, "VARIABLES= X, Y, UU UV UU_P UV_P \n");
    fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', rows, cols);
    dlmwrite(filename, RECSTATS, '-append', 'delimiter', ' ');
    
    disp('<strong>Phase-averaged velocity field exported</strong>');
    elseif tec==0
end

for i= 1:n_bins
    pcolor(reshape(phase_averaged_u{i},[rows,cols])'); shading interp; colormap jet; axis equal
    pause(0.5)
end

for i= 1:10
    pcolor(v(:,:,i)'); shading interp; colormap jet; axis equal
    pause(0.5)
end

tec = input('Write data to TECPLOT? [0/1] ');
if tec == 1
    disp('Working on results files...');
    fName = input('Enter file name to export: ', 's');
    fdName = input('Enter name of the folder: ', 's');
    mkdir(fdName);  % Create the folder if it doesn't exist

    for i = 1:n_bins
        stats = [x y phase_averaged_u{i} phase_averaged_v{i}];
        filename = fullfile(fdName, [fName, '_', num2str(i, '%04.f'), '.dat']);  % Full path to the file
        fid = fopen(filename, 'w');
        % Write header information
        fprintf(fid, 'TITLE=%s\n', filename);
        fprintf(fid, 'VARIABLES= x, y, U_phase, V_phase\n');
        fprintf(fid, 'ZONE I=%d J=%d F=POINT\n', rows, cols);
        disp(['Exporting data: ' num2str(i) ]);
        % Write data
        dlmwrite(filename, stats, '-append', 'delimiter', ' ');
        fclose(fid);
        
    end

    disp('DATA EXPORTED');
elseif tec == 0
    disp('Export cancelled.');
end 

    
    

