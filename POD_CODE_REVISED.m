clc; clear; close all;

pathToFolder = 'HCG050_NEAR_ANALYZED';
files = dir( fullfile(pathToFolder,'*.dat') );
NFiles = numel(files);
%read all files
dataV = cell(NFiles,1); %set up the length of the dataV

sumU= 0.0;
sumV= 0.0;
% Initialize the sum matrix
%% READ THE AND ARRANGE THEM IN X x Yx T 3D ARRAY 
disp('<strong>Reading files...</strong>')
for i=1:NFiles
    fid = fopen(fullfile(pathToFolder,files(i).name), 'rt');
    H = textscan(fid, '%s', 5, 'Delimiter','\n');  % Header lines to skip
    lineformat = repmat('%f',1,5); % Number of unique columns
    C = textscan(fid, lineformat, 'delimiter', ',','HeaderLines', 1, 'CollectOutput',1);
    fclose(fid);
    
    % Data columns for Uinst  and Vinst
    dataV{i} = C{1,1};
    x= dataV{i}(:,1);
    y= dataV{i}(:,2);
    yN= unique(C{1,1}(:,2));
    xN= unique(C{1,1}(:,1));
    rows= size(xN,1);
    cols= size(yN,1);
    dataV{i} = [C{1,1}(:,3) C{1,1}(:,4)];
    sumU= sumU + dataV{i}(:,1);
    sumV= sumV + dataV{i}(:,2);
    
end 
Umean= sumU/NFiles;
Vmean= sumV/NFiles;
disp('<strong>MEAN COMPUTED!</strong>')

% Remove the mean from the data
k= 0;
for j= 1:NFiles
    k= k+1;
    ufluc= dataV{k}(:,1) - Umean;
    vfluc= dataV{k}(:,2) - Vmean;
    
    UFLUC= reshape(ufluc,[rows,cols]);
    VFLUC= reshape(vfluc,[rows,cols]);

    U(:,:,k)= cat(3,UFLUC); %  U(:,i)= cat(3,ufluc) gives (rowsxcols,nFiles) % Don't need to do U_R reshape
    V(:,:,k)= cat(3,VFLUC); %  V(:,i)= cat(3,vfluc)
    
end

%% ANIMATE FEW FRAMES
% for ti=1:10
%     pcolor(X,Y,squeeze(U(:,:,ti)))
%     axis equal tight, shading interp
%     colormap('jet')
%     xlabel('X'), ylabel('Y')
%     pause(0.05)
%     drawnow
% end
disp('<strong>Reshaping Data...</strong>')
X= reshape(x,[rows,cols]);
Y= reshape(y,[rows,cols]);
U_R= reshape(U,rows*cols,size(U,3));
V_R= reshape(V,rows*cols,size(V,3));
vel_data= [U_R;V_R];
disp('<strong>COMPUTING MODAL COEFFICIENTS...</strong>')

%% POD ALGORTIHM
[psi,sigma,phi]= svd(U_R,'econ');
mode_energy= diag(sigma).^2;

%% ENERGY THRESHOLD
energy_threshold = 0.95; % Adjust as needed
cumulative_energy = cumsum(mode_energy) / sum(mode_energy);
num_modes_to_retain = find(cumulative_energy >= energy_threshold, 1);

POD_modes = psi; % (:, 1:num_modes_to_retain)
temporal_coeffs =  sigma*phi; %sigma(1:num_modes_to_retain, 1:num_modes_to_retain) * phi(:, 1:num_modes_to_retain)';

figure;
plot(cumulative_energy(1:10)*100, 'ko--');
xlabel('POD Modes');
ylabel('Cumulative Energy(%)');
% ylim([10,100])

percent_energy(1:10)= (mode_energy(1:10)/(sum(mode_energy)))*100;
figure(2);
plot(percent_energy(1:10),'square -')
xlabel('POD Modes'); ylabel('POD Energy (%)')
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
            contourf((reshape(POD_modes(:, i), [rows, cols])'), 20, 'LineColor', 'none');
            title(['POD Mode ', num2str(i)]);
            colorbar; colormap('jet')
        end
end
    


tec= input('Write data to TECPLOT? [0/1]');
if tec==1
    %% Export to TECPLOT format 
    disp('Working on results files...');
    PIVStats = [x y Umean Vmean POD_modes(:,1) POD_modes(:,2) POD_modes(:,3) POD_modes(:,4) POD_modes(:,5) POD_modes(:,6) POD_modes(:,7)];
    fName= input('Enter file name to export: ',"s");
    filename = [fName,'.dat'];
    fid = fopen(filename, 'w');
    fprintf(fid, 'TITLE=%s\n', filename);
    fprintf(fid, "VARIABLES= X, Y, U, V MODE1 MODE2 MODE3 MODE4 MODE5 MODE6 MODE7 \n");
    fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', rows, cols);
    dlmwrite(filename, PIVStats, '-append', 'delimiter', ' ');
    fclose(fid);
    
    disp('<strong>POD MODES EXPORTED!</strong>');
elseif tec==0
end

%% RECONSTRUCT THE DATA
r= 5; % reconstruct using r mode
U_red= psi(:,1:r) * sigma(1:r,1:r) * phi(:,1:r)';
t= 5; % no. instant
U_REDUCED= reshape(U_red(:,t),[rows,cols]);
contourf(U_REDUCED', 20, 'LineColor','none');
title(['Reconstructed Field using ',num2str(r),' modes'])
colormap('jet')
colorbar;
