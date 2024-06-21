% Developed by Sabal Bista (sblbista)
tic()
clear; clc; format long;
% Specify the directory where your CSV files are located
directory = './TEST2';

% List all CSV files in the directory
fileList = dir(fullfile(directory, '*.csv'));
u= cell(1,numel(fileList));
v= cell(1,numel(fileList));

for i = 1:numel(fileList)
    filename = fullfile(directory, fileList(i).name);
    
    % Read the CSV file, skipping the header row
    data = readmatrix(filename); % Skip 1 row (header) and 0 columns
    u{i}= data(:,2); v{i}= data(:,3);
    % Extract x, y, and z coordinates
    x = data(:, 11); % x coordinates 
    y = data(:, 12); % y coordinates 
    z = data(:, 13); % z coordinates 
    disp(['Loading timestep: ' num2str(i) ]);
end

vertices= cat(2,x,y);

% Select the range of region 
X_min = 0.0;
X_max = 0.254;
Y_min = 0;
Y_max = 0.1016;

% Triangulation
DT= delaunay(vertices);

u_all= cell2mat(u); 
v_all= cell2mat(v);
% Plot the triangulated data

% figure(2)
% trisurf(DT,vertices(:,1),vertices(:,2),u_all(:,1),'EdgeColor','none'); 
% axis tight; daspect([1 1 1]); 
% shading interp; colormap gray; view(2);title('Velocity U(i)'); colorbar;
% figure(3)
% trisurf(DT,vertices(:,1),vertices(:,2),v_all(:,1),'EdgeColor','none'); 
% axis tight; daspect([1 1 1]); 
% shading interp; colormap gray; view(2);title('Velocity U(j)'); colorbar;

% Save MAT file
save_inst= input('Save instantenous data in a MAT file? [0/1] ');
if save_inst==1
    save('inst_dataCFD.mat','u_all','v_all','x','y','-v7.3');
    data= load('inst_dataCFD.mat');
elseif save_inst==0
end

disp('<strong>INTERPOLATING DATA...</strong>')
%% Interpolate data in a cartesian grid
x_round=(round(x,4)); % round x and y to 4 digits
y_round=(round(y,4)); %
x_round= x_round(x_round>=X_min & x_round<=X_max);
y_round= y_round(y_round>=Y_min & y_round<=Y_max);

I= 250; J= 100; % assign no. of grids for interpolation 
x1= linspace(min(x_round),max(x_round),I);
y1= linspace(min(y_round),max(y_round),J);

[xq,yq]= meshgrid(x1,y1);

[nx ny]=size(xq);
F= cell(1,numel(fileList));


for k=1:numel(fileList)

    F1{k}= scatteredInterpolant(vertices(:,1),vertices(:,2),u_all(:,k)); % griddata(vertices(:,1),vertices(:,2),u_all(:,k),xq,yq);
    u_grid= F1{k}(xq,yq);
    u_temp1(:,k)= reshape(u_grid,[I*J,1]);
    u_temp(:,:,k)= cat(3,u_grid);

    F2{k}= scatteredInterpolant(vertices(:,1),vertices(:,2),v_all(:,k));
    v_grid= F2{k}(xq,yq);
    v_temp(:,:,k)= cat(3,v_grid);
    v_temp1(:,k)= reshape(v_grid,[I*J,1]);
    
    disp(['Interpolating timestep: ' num2str(k) ]);
    
end

X= reshape(xq,[nx*ny,1]);
Y= reshape(yq,[nx*ny,1]);

%% ANIMATE FEW FRAMES
animate= input('Animate few frames? [0/1] ');
if animate== 1
    for ti=1:10
        pcolor(xq,yq,squeeze(u_temp(:,:,ti)))
        axis equal tight, shading interp
        colormap('gray')
        xlabel('X'), ylabel('Y')
        pause(0.1)
        drawnow
    end
elseif animate==0
end

% Save interpolated MAT file
save_inst= 1; % input('Save interpolated instantenous data in a MAT file? [0/1] ');
if save_inst==1
    save('inst_data_interpolated.mat','u_temp','v_temp','xq','yq','-v7.3');
    % data= load('inst_data_interpolated.mat');
elseif save_inst==0
end

% %% Calculate mean
sumUi= 0; sumUj= 0; nFiles= size(u_temp,3);
for i= 1:nFiles
    sumUi= sumUi + u_temp1(:,i);
    sumUj= sumUj + v_temp1(:,i);
    i=i+1;
end
Umean= sumUi/nFiles;
Vmean= sumUj/nFiles;

UMEAN= reshape(Umean,nx,ny);
figure(2)
% contourf(xq,yq,UMEAN,40,'LineColor','none'); 
pcolor(xq,yq,UMEAN); colormap gray; axis equal tight, shading interp

%% Calculate the fluctuating components
for j= 1:nFiles
    uFluc(:,j)= u_temp1(:,j) - Umean;
    vFluc(:,j)= v_temp1(:,j) - Vmean;
    disp(['Calculating the fluctuating component of timestep: ' num2str(j) ]);
end

%% Calculate the modes 
% ---------------------------------------------------------------------
disp('Performing SVD...')
[psi,sigma,phi]= svd(uFluc,'econ');
mode_energy= diag(sigma).^2;
POD_modes = psi;

percent_energy(1:10)= (mode_energy(1:10)/(sum(mode_energy)))*100;
figure(2);
plot(percent_energy(1:10),'square -')
xlabel('POD Modes'); ylabel('POD Energy')
title('Energy distributions among the first 10 POD modes')

visMod= input('Visualize modes? [0/1]');
if visMod== 1
    nModes= input('Enter number of modes to visualize(>3): ');
    p= ceil(nModes/2);
    q= ceil(nModes/3);

        for i = 1:nModes
            figure(2);
            % hold(q,p,i),'on');
            subplot(q,p,i); 
            contourf((reshape(POD_modes(:, i), [nx, ny])), 20, 'LineColor', 'none');
            title(['POD Mode ', num2str(i)]);
            colorbar; colormap('jet')
            axis('equal')
        end
end

%% Export to TECPLOT format 
tec= input('Write data to TECPLOT? [0/1]');
if tec==1
    disp('Working on results files...');
    PIVStats = [X Y Umean Vmean psi(:,1) psi(:,2) psi(:,3) psi(:,4) psi(:,5) psi(:,6) psi(:,7) psi(:,8) psi(:,9) psi(:,10)];
    fName= input('Enter file name to export: ',"s");
    filename = [fName,'.dat'];
    fid = fopen(filename, 'w');
    fprintf(fid, 'TITLE=%s\n', filename);
    fprintf(fid, "VARIABLES= X, Y, U, V, MODE1, MODE2, MODE3, MODE4 MODE5, MODE6, MODE7, MODE8 MODE9 MODE10\n");
    fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', nx, ny);
    dlmwrite(filename, PIVStats, '-append', 'delimiter', ' ');
    fclose(fid);

    disp('<strong>POD MODES EXPORTED!</strong>');
elseif tec==0 
end
toc()
