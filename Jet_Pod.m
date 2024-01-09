%
%% Interpolates CFD data in a cartesian grid and performs SVD to get POD information (Also writes the data for the SPOD code) %%
%
clear,clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathToFolder = 'Instant_Quant';
files = dir( fullfile(pathToFolder,'*.csv') );
%# read all files
dataV = cell(numel(files),1); %set up the length of the dataV
count= 0;
sumUi= 0.0; 
sumUj= 0.0;
sumUk= 0.0;
sumVf= 0.0;
sumPr= 0.0; % Pressure
nfile= numel(files);
skip= 1; % no. of files to skip (set to 1 if no skip)
i=0; 
for k=1:skip:nfile
    i= i+1;
    % time(i)= (i)/f;
    count= count+1;
    fid= fopen(fullfile(pathToFolder,files(i).name), 'rt');
    H= textscan(fid, '%s', 1, 'Delimiter','\n');
    %H store the data and ignore the first header line
    %C = textscan(fid, repmat('%f ',1),repmat('%f ',1), 'Delimiter',' ', ...
    %'MultipleDelimsAsOne',true, 'CollectOutput',true);    
    lineformat = repmat('%f',1,14);
    C1 = textscan(fid, lineformat, 'delimiter', ',','HeaderLines', 0, 'CollectOutput',1);
    fclose(fid);
    %XY plane replace nans with zeros
    C = sortrows(C1{1,1},[12,13]);
    [m,n]=size(C);
    k = (isnan(C));
    for ii=1:m
        for jj=1:n
            if (k(ii,jj)==1)
               C(ii ,jj)=0;
            end
        end
    end	
    dataV{i} = [C(:,3) C(:,4) C(:,5) C(:,12) C(:,13) C(:,14)]; % u v w x y z

    sumUi=sumUi+dataV{i}(:,1);
    sumUj=sumUj+dataV{i}(:,2);
    sumUk=sumUk+dataV{i}(:,3);
    %digits(4)
    x= dataV{i}(:,5);
    y= dataV{i}(:,6);

end

% Calculate mean quantities
meanUi=sumUi/count;
meanUj=sumUj/count;
meanUk=sumUk/count;
disp('<strong>MEAN CALCULATED!</strong>');

%% Grid Data
x_temp=(round(C(:,12),4)); % round x and y to 4 digits
y_temp=(round(C(:,13),4)); %

x_grid_no= 100; y_grid_no= 100; % assign no. of grids for interpolation 
x1= linspace(min(x_temp),max(x_temp),x_grid_no);
y1= linspace(min(y_temp),max(y_temp),y_grid_no);

[xq,yq]= meshgrid(x1,y1);

[nx ny]=size(xq);

disp('Calculating fluctuating components')

k=0;
for i=1:skip:nfile
    k=k+1;
    ud{k} = dataV{k}(:,1)-meanUi;
    vd{k} = dataV{k}(:,2)-meanUj;  
    % wd = dataV{k}(:,3)-meanUk;
    u_grid= griddata(C(:,12),C(:,13),ud{k}(:,1),xq,yq); % interpolates u and v to new grid
    v_grid= griddata(C(:,12),C(:,13),vd{k}(:,1),xq,yq);
    ud_temp(:,k)= reshape(u_grid,[nx*ny,1]); %cat(3,u_grid);
    vd_temp(:,k)= reshape(v_grid,[nx*ny,1]);
    veldata(:,k)=[ud_temp(:,k); vd_temp(:,k)];
end
X= reshape(xq,[nx*ny,1]);
Y= reshape(yq,[nx*ny,1]);

%% Interpolated mean
UMean_int= griddata(C(:,12),C(:,13),meanUi,xq,yq);
VMean_int= griddata(C(:,12), C(:,13),meanUj,xq,yq);
UMean_int= reshape(UMean_int,[nx*ny,1]);
VMean_int= reshape(VMean_int,[nx*ny,1]);

cord=[X Y];

% ---------------------------------------------------------------------
%% Calculate the modes 
% ---------------------------------------------------------------------
disp('Performing SVD...')
[psi,sigma,phi]= svd(vd_temp,'econ');
mode_energy= diag(sigma).^2;
POD_modes = psi;

percent_energy(1:10)= (mode_energy(1:10)/(sum(mode_energy)))*100;
figure(2);
plot(percent_energy(1:10),'square -')
xlabel('POD Modes'); ylabel('POD Energy (%)')
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
disp('Working on results files...');
PIVStats = [X Y UMean_int VMean_int psi(:,1) psi(:,2) psi(:,3) psi(:,4) psi(:,5) psi(:,6) psi(:,7) psi(:,8) psi(:,9) psi(:,10)];
fName= input('Enter file name to export: ',"s");
filename = [fName,'.dat'];
fid = fopen(filename, 'w');
fprintf(fid, 'TITLE=%s\n', filename);
fprintf(fid, "VARIABLES= X, Y, U, V, MODE1, MODE2, MODE3, MODE4 MODE5, MODE6, MODE7, MODE8 MODE9 MODE10\n");
fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', nx, ny);
dlmwrite(filename, PIVStats, '-append', 'delimiter', ' ');
fclose(fid);

disp('<strong>POD MODES EXPORTED!</strong>');

%% Export data to SPOD readable format (nt,nx,ny)
for l= 1:length(ud_temp(1,:))
    uSPOD(:,:,l)= reshape(ud_temp(:,l),[nx,ny]);
    vSPOD(:,:,l)= reshape(vd_temp(:,l),[nx,ny]);
end
uSPOD= permute(uSPOD,[3,1,2]); vSPOD= permute(vSPOD,[3,1,2]); 

save('v_velocity.mat','uSPOD','vSPOD','xq','yq', '-v7.3')
% save cord.mat X Y
