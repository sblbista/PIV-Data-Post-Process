% This code reads PIV data and computes the POD modes, energy distrubution,
% and mode coefficients using SVD technique. It also reconstructs the
% velocity field using 50% energy criteria and displays the velocity
% magnitude. It also computes the phase angles based on the first two
% dominant modes. 

clc; clear; close all;

pathToFolder = 'HCG200_NWAKE_PROCESSED';
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
    disp(['Time instance read: ' num2str(i) ]);
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
grid_points= rows*cols;
U_R= reshape(U,rows*cols,size(U,3));
V_R= reshape(V,rows*cols,size(V,3));
vel_data= [U_R;V_R];
disp('<strong>COMPUTING MODAL COEFFICIENTS...</strong>')

%% POD ALGORTIHM
[psi,sigma,phi]= svd(vel_data,'econ');
mode_energy= diag(sigma).^2;

%% ENERGY THRESHOLD
energy_threshold = 0.50; % Adjust as needed
cumulative_energy = cumsum(mode_energy) / sum(mode_energy);
num_modes_to_retain = find(cumulative_energy >= energy_threshold, 1);

POD_modes = psi; % (:, 1:num_modes_to_retain)
temporal_coeffs =  sigma*phi; %sigma(1:num_modes_to_retain, 1:num_modes_to_retain) * phi(:, 1:num_modes_to_retain)';

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
            contourf((reshape(POD_MAG(1:[rows*cols], i), [rows, cols])'), 20, 'LineColor', 'none');
            title(['POD Mode ', num2str(i)]);
            colorbar; colormap('jet')
        end
end
   
A= vel_data'*POD_modes;

tec= input('Write data to TECPLOT? [0/1]');
if tec==1
    %% Export to TECPLOT format 
    disp('Working on results files...');
    PIVStats = [x y Umean Vmean PODU(:,1) PODU(:,2) PODU(:,3) PODU(:,4) PODU(:,5) PODV(:,1) PODV(:,2) PODV(:,3) PODV(:,4) PODV(:,5)];
    fName= input('Enter file name to export: ',"s");
    filename = [fName,'.dat'];
    fid = fopen(filename, 'w');
    fprintf(fid, 'TITLE=%s\n', filename);
    fprintf(fid, "VARIABLES= X, Y, U, V MODE1U MODE2U MODE3U MODE4U MODE5U MODE1V MODE2V MODE3V MODE4V MODE5V \n");
    fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', rows, cols);
    dlmwrite(filename, PIVStats, '-append', 'delimiter', ' ');
    fclose(fid);
    
    disp('<strong>POD MODES EXPORTED!</strong>');
elseif tec==0
end

%% RECONSTRUCT THE DATA
r= num_modes_to_retain; % reconstruct using r mode
U_red= PODU(:,1:r) * sigma(1:r,1:r) * phi(:,1:r)'; %psi(:,1:r) * sigma(1:r,1:r) * phi(:,1:r)';
V_red= PODV(:,1:r) * sigma(1:r,1:r) * phi(:,1:r)'; 
t= 1; % no. instant
U_REDUCED= reshape(Umean,[rows,cols]) + reshape(U_red(:,t),[rows,cols]);
V_REDUCED= reshape(Vmean,[rows,cols]) + reshape(V_red(:,t),[rows,cols]);
UV_REDUCED_MAG= sqrt(U_REDUCED.^2 + V_REDUCED.^2);
contourf(UV_REDUCED_MAG', 20, 'LineColor','none');
title(['Reconstructed Field using ',num2str(r),' modes'])
colormap('jet')
colorbar;
axis equal

% Plot the time coefficients
% a1= phi(:,1); a2= phi(:,2);
% a1= temporal_coeffs(:,1); a2= temporal_coeffs(:,2);
a1= A(:,1); a2= A(:,2);
lambda1= mode_energy(1)/(sum(mode_energy)); lambda2= mode_energy(2)/(sum(mode_energy));
plot(a1,a2,'xk')
axis equal
xlabel('a_1'), ylabel('a_2')
for i= 1:length(a1)
    r(i)= sqrt(a1(i)^2 + a2(i)^2);
    a1S(i)= a1(i)/r(i); 
    a2S(i)= a2(i)/r(i); 
end
plot(a1S,a2S,'o'); axis equal; xlabel('a_1*');ylabel('a_2*')

plot(a1/sqrt(2*lambda1), a2/sqrt(2*lambda2), 'ro')
axis equal
xlabel('$a_1(t)/\sqrt{2 \lambda_1}$','Interpreter','latex')
ylabel('$a_2(t)/\sqrt{2 \lambda_2}$','Interpreter','latex')

phase_angle= rad2deg(atan2(a1*sqrt(2*lambda2), a2*sqrt(2*lambda1)));

phase_angle= sort(phase_angle,'ascend');

