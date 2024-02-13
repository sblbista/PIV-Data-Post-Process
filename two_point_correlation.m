%% Spatial Two-point correlation
clc; clear; close all;

pathToFolder = 'HCG050_NEAR_ANALYZED'; % HCG050_NEAR_ANALYZED
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

xref= 0.0198694276338464842; yref= 0.0453171724133631243; % Get from Tecplot

[d, x_idx]= min(abs(xN-xref)); 
[d, y_idx]= min(abs(yN-yref));

% Remove the mean from the data
k= 0;

UPrime2= 0; 
VPrime2= 0;
Ruu1= 0; RUU2= 0;
Rvv1= 0; RVV2= 0;
for j= 1:NFiles
    k= k+1;
    ufluc= dataV{k}(:,1) - Umean;
    vfluc= dataV{k}(:,2) - Vmean;
    
    UFLUC= reshape(ufluc,[rows,cols]);
    VFLUC= reshape(vfluc,[rows,cols]);

    U(:,:,k)= cat(3,UFLUC); %  U(:,i)= cat(3,ufluc) gives (rowsxcols,nFiles) % Don't need to do U_R reshape
    V(:,:,k)= cat(3,VFLUC); %  V(:,i)= cat(3,vfluc)

    UPrime2 =(ufluc).^2 + UPrime2;
    VPrime2 =(vfluc).^2 + VPrime2; 
    urms = sqrt(UPrime2/NFiles);
    vrms = sqrt(VPrime2/NFiles);
    URMS= reshape(urms,[rows,cols]);
    VRMS= reshape(vrms,[rows,cols]);

    % Two-point spatial Correlation
    Ruu1= Ruu1 + (UFLUC(x_idx,y_idx).*UFLUC);
    RUU1= (Ruu1)/NFiles;
    RUU2= (URMS(x_idx,y_idx).*URMS);
    Ruu= RUU1./RUU2;
    
    Rvv1= Rvv1 + (VFLUC(x_idx,y_idx).*VFLUC);
    RVV1= (Rvv1)/NFiles;
    RVV2= (VRMS(x_idx,y_idx).*VRMS);
    Rvv= RVV1./RVV2;
end



Ruu= reshape(Ruu,[rows*cols,1]);
Rvv= reshape(Rvv,[rows*cols,1]);
toremove= isnan(Ruu);
Ruu(toremove)= 0;
Rvv(toremove)= 0;
disp('<strong>RMS CALCULATION COMPLETE!</strong>');
disp('<strong>SPATIAL TWO-POINT CORRELATION CALCULATION COMPLETE!</strong>');


tec= input('Write data to TECPLOT? [0/1]');
if tec==1
    %% Export to TECPLOT format 
    disp('Working on results files...');
    PIVStats = [x y Umean Vmean Ruu Rvv urms vrms];
    fName= input('Enter file name to export: ',"s");
    filename = [fName,'.dat'];
    fid = fopen(filename, 'w');
    fprintf(fid, 'TITLE=%s\n', filename);
    fprintf(fid, "VARIABLES= X, Y, U, V R<sub>uu</sub> R<sub>vv</sub> URMS VRMS \n");
    fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', rows, cols);
    dlmwrite(filename, PIVStats, '-append', 'delimiter', ' ');
    fclose(fid);
    
    disp('<strong>POD MODES EXPORTED!</strong>');
elseif tec==0
end

disp('<strong>CODE ENDED!</strong>');

