clear; clc
% file_name= 'Large_field_XY_3mps_stats.nc';
% 
% ncdisp(file_name);
% grid_y = ncread(file_name,'grid_y');
% grid_x = ncread(file_name,'grid_x');
% 
% i=10;
% vel_x_mean = ncread(file_name,'vel_x_mean');
% vel_y_mean = ncread(file_name,'vel_y_mean');
% 
% vel_x_rms = ncread(file_name,'vel_x_rms');
% vel_y_rms = ncread(file_name,'vel_y_rms');
% Re_xy = ncread(file_name,'Re_xy');
M= load("data1.dat");
y= M(:,1);
UUe= M(:,2);

% flag= input('Write data to TECPLOT? [0/1]');
% 
% if flag== 1
% %% Exporting to TECPLOT format
% 
%     Umean= reshape(vel_x_mean,1,[])';
%     Vmean= reshape(vel_y_mean,1,[])';
%     X= reshape(grid_x,1,[])';
%     Y= reshape(grid_y,1,[])';
%     
%     % Results 
%     disp('Working on results files...');
%     PIVStats = [X Y Umean Vmean];
%     filename = 'MPSNew.dat';
%     fid = fopen(filename, 'w');
%     fprintf(fid, 'TITLE=%s\n', filename);
%     fprintf(fid, 'VARIABLES= X, Y, U, V\n');
%     fprintf(fid, 'ZONE  I=%d  J=%d F=POINT\n', 766, 199);
%     dlmwrite(filename, PIVStats, '-append', 'delimiter', ' ');
%     fclose(fid);
% elseif flag== 0
% end

% % Write velocity profile 
% profileHigh1= [vel_x_mean(700,:)' grid_y(1,:)'];
% csvwrite('profile@High1.csv',profileHigh1);

%% Double Averaging
% y= grid_y(1,:)';
% UDoubleAvg= mean(vel_x_mean,1)';
% VDoubleAvg= mean(vel_y_mean,1)';
% UVDoubleAvg= mean(Re_xy,1)';
% uu= vel_x_rms.^2;
% UUDoubleAvg= mean(uu,1)';
% gradUY= gradient(UDoubleAvg, y);
% Initial guess for uTau
nu= 1.5e-5;
% tau= (nu*gradUY - UVDoubleAvg);


% plot(tau,y)
% xlabel('u_{\tau}')
% ylabel('y')
% 
% plot(-UVDoubleAvg,y)
% xlabel("-u'v'")
% ylabel('y')

%% Extrapolate utau
% Linear fit of Reynolds stress
% p= polyfit(-UVDoubleAvg,y,1);
% 
% utau1= sqrt(-p(2)/p(1));
% fprintf('u* from Reynolds stress= %f m/s \n',utau1)
% 
% % Linear fit of total Shear Stress
% p2= polyfit(tau,y,1);
% 
% utau2= sqrt(-p2(2)/p2(1));
% fprintf('u* from total shear stress= %f m/s \n',utau2)
% 
% %% Compare with analytical plot
% k= 0.41; B= 5.0;
% vv= linspace(5,1000);
% y11= (1/k) * log(vv) + B;
% 
% semilogx(vv,y11,'-k')
% xlabel('y^+')
% ylabel('U^+')
% hold on 
% 
% % Log law plot from Reynolds shear stress
% yP1= (y*utau1)/nu;
% uP1= UDoubleAvg/utau1;
% 
% semilogx(yP1,uP1,'ob');
% 
% % Log law plot from Total Shear stress
% yP2= (y*utau2)/nu;
% uP2= UDoubleAvg/utau2;
% 
% semilogx(yP2,uP2,'square','MarkerFaceColor','r')
% hold off
% 
% %% New method
% T= (nu*gradUY - UVDoubleAvg);
% U= 3; % m/s
% % Calculate the displacement thickness delta*
% dispThickness= trapz(y,(1 - (UDoubleAvg/U)));
% 
% % Calculate the momentum thickness theta
% momentumThickness= trapz(y,(UDoubleAvg/U).*(1 - (UDoubleAvg/U)));
% 
% % Calculate the shape factor (H)
% H= dispThickness/momentumThickness;
% 
% eta= 0:0.01:1;
% a= 0.5055; b= 1.156;
% VVe= VDoubleAvg/max(VDoubleAvg);%tanh(a*eta + b*eta.^3);
% plot(y,VVe)
% utau3= sqrt(T(1)./ (H .* (1-VVe(1)) + (H-1).*(eta-1)));

% plot(eta,utau3,'o');


%% Optimizing Krogstad's equation

delta= 0.01835; %0.3; % Boundary layer thickness 30 cm
yDelta= y./delta;

% fPI= @(PI) (U - UDoubleAvg)./(utau2) - (2.*PI./k) .* (1 - (0.5.*PI) .* ((1 + 6.*PI) .* (yDelta).^2 - (1 + 4.*PI).*(yDelta).^3)) - (1./k).*log(yDelta);

k= 0.41;
Ue= 33.32;
u= Ue*UUe;
A= (-2 + 6 * yDelta.^2 - 4 * yDelta.^3);
B= (yDelta.^2 - yDelta.^3 - log(yDelta) - k.*((Ue - u)/1.5));

fPI= @(PI) A.*PI + B;

PIopt= - (A.'*B) / (A.'*B) 



