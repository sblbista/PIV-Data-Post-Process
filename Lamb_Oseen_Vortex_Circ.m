% Lamb-Oseen vortex
clc; clear; close all;
r=-3.1:0.15:3.1;

[x,y]=meshgrid(r);
dx= 0.15; dy= 0.15;
rr= x.^2+y.^2;
u= -y./rr.*(1-exp(-rr));
v= x./rr.*(1-exp(-rr));
quiver(x,y, u,v , "b" ) ;
hold on;
axis([-3 3 -3 3]);

axis equal ;

%% Calcualte Velocity Gradient
[dudx, dudy]= gradient(u,dx,dy);
[dvdx, dvdy]= gradient(v,dx,dy);

wz= (dvdx-dudy);
%contour(x,y,wz)

for i= 1:numel(x(:,1))
    for j= 1:numel(y(:,1))
        D2D= [dudx(i,j) dudy(i,j); dvdx(i,j) dvdy(i,j)];
        temp_eig= eig(D2D);
        lambda(i,j)= imag(temp_eig(1,1));
        lambda1(i,j)=lambda(i,j)*sign(wz(i,j));
        P(i,j)= trace(D2D);
        Q(i,j)= -(dudy(i,j).*dvdx(i,j)) + (dudx(i,j).*dvdy(i,j));
    end
end

contour(x,y,lambda1)

X= reshape(x,[numel(x),1]);
Y= reshape(y,[numel(x),1]);
U= reshape(u,[numel(x),1]);
V= reshape(v,[numel(x),1]);
Wz= reshape(wz,[numel(x),1]);
ld1= reshape(lambda1,[numel(x),1]);
P= reshape(P,[numel(x),1]);
Q= reshape(Q,[numel(x),1]);
%% Export to TECPLOT format 
disp('Working on results files...');
PIVStats = [X Y U V Wz ld1 P Q];
filename = 'lambos.dat';
fid = fopen(filename, 'w');
fprintf(fid, 'TITLE=%s\n', filename);
fprintf(fid, "VARIABLES= X, Y, U, V, VORT LAMBDA1 P1 Q1\n");
fprintf(fid, 'ZONE  I= %d  J= %d F=POINT\n', 42, 42);
dlmwrite(filename, PIVStats, '-append', 'delimiter', ' ');
fclose(fid);