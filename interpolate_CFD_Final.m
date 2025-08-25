% Developed by Sabal Bista (sblbista)
tic()
clear; clc; format long;
% Specify the directory where your CSV files are located
directory = 'D:/Data/nasif/bistas/thakurta/IDDES/PC_data/Horizontal_midplane_inst/';


% List all CSV files in the directory
fileList = dir(fullfile(directory, '*.csv'));
skip= 1;
nFiles= floor(numel(fileList)./skip); 
 
u= []; v= [];
u_grid= []; v_grid= []; 
u_temp= []; v_temp= []; 
u_column= 2; v_column= 4; 
x_column= 10; y_column= 12; 

d= 0.0255; % length scale of the body
grid_size= 0.001; % 1mm
% Load data
for i = 1:50
    filename = fullfile(directory, fileList(i).name);
    
    % Read the CSV file, skipping the header row
    data = readmatrix(filename); % Skip 1 row (header) and 0 columns

    % Remove rows with NaN values in any column (except for x and y)
    validData = all(~isnan(data(:, [u_column, v_column])), 2); % Check if u, v, w, p are not NaN
    data = data(validData, :); % Keep only valid rows

    u= data(:,u_column); v= data(:,v_column);
    % x and y coordinates
    x = data(:,x_column); % x coordinates 
    y = data(:,y_column); % y coordinates 
    
    % Select the range of region 
    X_min = -1*d;
    X_max = 8*d;
    Y_min = -2*d;
    Y_max = 2*d;
    
    x_round=(round(x,4)); % round x and y to 4 digits
    y_round=(round(y,4)); %
    x_round= x_round(x_round>=X_min & x_round<=X_max);
    y_round= y_round(y_round>=Y_min & y_round<=Y_max);
    
    I= round((abs(X_max) + abs(X_min))/(grid_size)); J= round((abs(Y_max) + abs(Y_min))/(grid_size)); % assign no. of grids for interpolation 
    x1= linspace(min(x_round),max(x_round),I);
    y1= linspace(min(y_round),max(y_round),J);
    
    [xq,yq]= ndgrid(x1,y1);
    
    [nx, ny]=size(xq);

    % Interpolating data
    F1 = scatteredInterpolant(x,y,u,"natural"); 
    u_grid = F1(xq,yq);
    u_temp(:,i) = reshape(u_grid,[I*J,1]);
    
    F2 = scatteredInterpolant(x,y,v,"natural"); 
    v_grid = F2(xq,yq);
    v_temp(:,i) = reshape(v_grid,[I*J,1]);

        
    disp(['Loading timestep and interpolating: ' num2str(i) ]);
end


% % Save MAT file
% save_inst= input('Save instantenous data in a MAT file? [0/1] ');
% if save_inst==1
%     save('inst_dataCFD.mat','u_all','v_all','x','y','-v7.3');
%     data= load('inst_dataCFD.mat');
% elseif save_inst==0
% end


X= reshape(xq,[nx*ny,1]);
Y= reshape(yq,[nx*ny,1]);

%% ANIMATE FEW FRAMES
animate= 1; % input('Animate few frames? [0/1] ');
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
save_intr= input('Save interpolated instantenous data in a MAT file? [0/1] ');
if save_intr==1
    save('inst_data_interpolated.mat','u_temp','v_temp','w_temp','p_temp','xq','yq','-v7.3');
    % data= load('inst_data_interpolated.mat');
elseif save_inst==0
end
