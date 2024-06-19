%% 
% This code loads the CSV file exported from STAR-CCM+ and 
% triangulates the data using Delaunay method
%%
clear; clc;
% Specify the directory where your CSV files are located
directory = './TEST';

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
end

vertices= cat(2,x,y);

% Triangulation
DT= delaunay(vertices);


u_all= cell2mat(u); 
v_all= cell2mat(v);
figure(1)
trisurf(DT,vertices(:,1),vertices(:,2),u_all(:,1),'EdgeColor','none'); 
axis tight; daspect([1 1 1]); 
shading interp; colormap gray; view(2);title('Velocity U(i)'); colorbar;
figure(2)
trisurf(DT,vertices(:,1),vertices(:,2),v_all(:,1),'EdgeColor','none'); 
axis tight; daspect([1 1 1]); 
shading interp; colormap gray; view(2);title('Velocity U(j)'); colorbar;

