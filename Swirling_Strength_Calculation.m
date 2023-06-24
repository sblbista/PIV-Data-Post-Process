clc; clear; close all;


tic()
data= readtable('SCC000Z00B_Processed/Right_Processed_1950.dat'); % SCC000Z00B_Processed/Right_Processed_0117.dat Jet_0001.txt
data= data{:,1:6};
X= data(:,1);
Y= data(:,2);

U= data(:,3);
V= data(:,4);

x= unique(X);
y= unique(Y);

shapeX= size(x,1);
shapeY= size(y,1);

u=reshape(U,[shapeX,shapeY])';
v=reshape(V,[shapeX,shapeY])';
% Plot Contour of velocities
    
contourf(x,y,v,'LineColor','none'), colorbar; axis equal; 
figure(2); contourf(x,y,v,'LineColor','none'), colorbar;axis equal

dvdx= zeros(size(u));
dudy= zeros(size(u));
dvdy= zeros(size(u));
dudx= zeros(size(u));
dx= x(2) - x(1);
dy= y(2) - y(1);


[dudx, dudy]= gradient(u,dx,dy);
[dvdx, dvdy]= gradient(v,dx,dy);
wz= (dvdx - dudy);
% contourf(x,y,wz,8,'LineColor','none'), colorbar;
for i= 1:shapeY
    for j= 1:shapeX
        D2D= [dudx(i,j) dudy(i,j); dvdx(i,j) dvdy(i,j)];
        temp_eig= eig(D2D);
        lambda(i,j)= imag(temp_eig(1,1));
        lambda1(i,j)=lambda(i,j)*sign(wz(i,j));
        P(i,j)= -trace(D2D);
        Q(i,j)= -(dudy(i,j).*dvdx(i,j)) + (dudx(i,j).*dvdy(i,j));
 
    end
end

% Calculate RMS of LambdaCI
lambdaCI_SUM= 0;

for i= 1:size(U,1)
    lambdaCI_SUM= lambda1(i).^2 + lambdaCI_SUM;
end
lambdaCI_RMS= sqrt((1/size(U,1) * lambdaCI_SUM));
contourf(x,y,lambda1,'LineColor','none')

wz= wz';
lambda1= lambda1';

P= P'; Q= Q';
Wz= reshape(wz,size(U));
ld1= reshape(lambda1,size(U));
P= reshape(P,size(U));
Q= reshape(Q,size(U));

% Remove noise from data 
for i= 1:size(U,1)
    if abs(ld1(i))/lambdaCI_RMS <= 1.5
        ld1(i)= 0;
    end
end

flag= input('Write data as CSV file? [0/1]');

% Galilean Decomposition
% uG= u - 0.8* 0.312;
% uG= uG';
% UG= reshape(uG,size(U));

if flag==1
    %% Export to TECPLOT format 
    disp('Working on results files...');
    PIVStats = [X Y U V Wz ld1 P Q];
    filename = 'SCC_1950_LD1.dat';
    fid = fopen(filename, 'w');
    fprintf(fid, 'TITLE=%s\n', filename);
    fprintf(fid, "VARIABLES= X, Y, U, V, VORT LAMBDA1 P Q \n");
    fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', shapeX(1), shapeY(1));
    dlmwrite(filename, PIVStats, '-append', 'delimiter', ' ');
    fclose(fid);
end

toc()
