clc; clear; close all;
format long
M= load('john_FOVB_VelocityProfile.csv');


y= M(:,1);
%y= y/1000;
u= M(:,2);
dUdY= gradient(u,y);
% Plot the data
plot(y, u,'ok')
xlabel('y(m)')
ylabel('u(m/s)')


%% Plot log law with constant kapa=0.41 

% Change the y to log scale
flag= input('Take all the data for calculation? [y/n]','s');

if flag=='y'
    pomx= log(y);
    plot(pomx, u,'square')
    xlabel('log(y)');
    ylabel('U')
    legend('Data Points')

    % Perform least square fitting
    %p= polyfit(pomx(2:end),u(2:end),1);
    p= polyfit(pomx,u,1);
    % Intercepts
    a= p(1);
    b= p(2);
    hold on
    % Plot the fitted curve on the data
    plot(pomx,polyval(p,pomx))
    legend('Data Points','Fitted Curve')

elseif flag=='n'
    r1= input('Enter lower bound of data: ');
    r2= input('Enter upper bound of data: ');
    pomx= log(y(r1:r2));
    u1= u(r1:r2);
    plot(pomx, u1,'square')
    xlabel('log(y)');
    ylabel('U')
    legend('Data Points')

    % Perform least square fitting
    p= polyfit(pomx,u1,1);
    % Intercepts
    a= p(1);
    b= p(2);
    hold on
    % Plot the fitted curve on the data
    plot(pomx,polyval(p,pomx))
    legend('Data Points','Fitted Curve')

end

hold off
%% Calculate the friction velocity and B from log law

kappa= 0.41; % von Karman constant
nu= 1.0e-6; % kinematic viscosity m2/s

% Friction velocity
uStar= a*kappa

% log law constant
B= (b./uStar)-((1/kappa).*log(uStar/nu))

%% Plot log law with the calcualted values

ypl= (y.*uStar)./nu;
Upl= (u./uStar);
semilogx(ypl,Upl,'r*')
xlabel('y+')
ylabel('U+')

hold on

% Comapre with analytical log law
vv=linspace(5,1000);
yll=1/kappa.*log(vv)+5.0;
semilogx(vv,yll,'b')
xlabel('y+')
ylabel('U+')

legend('log law fit','Analytical Log law')

