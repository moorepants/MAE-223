% File: bicycle.m
% Date: July 17, 2006
% Author: Jason Moore
% Description:
clear
close all
%--------Vary Front Wheel Diameter
lambda  = 73*pi/180;                     % head tube angle [rad]
tr      = 0.0602;                        % trail [m]
wb      = 1.0287;                        % wheel base [m]
dfmin = 0.32;
dfmax = 0.78;
dfint = 0.02;
n = ((dfmax-dfmin)/dfint)+1;
for i = 1:n
    df(i) = (i-1)*dfint+dfmin;
    vdata(i,1:2) = v_crit(df(i),lambda,tr,wb);
end
figure(1)
plot(df,vdata(:,1),'r',df,vdata(:,2),'b')
xlabel('Front Wheel Diameter [m]')
ylabel('Velocity [m/s]')
legend('Weave Critical Velocity','Capsize Critical Velocity')
%---------Vary Head Tube Angle
clear
df      = 0.6985;                        % front wheel diameter [m]
tr      = 0.0602;                        % trail [m]
wb      = 1.0287;                        % wheel base [m]
lambdamin = 1.10; % 63*pi/180;
lambdamax = 1.48; % 85*pi/180;
lambdaint = 0.02;
n = ((lambdamax-lambdamin)/lambdaint)+1;
for i = 1:n
    lambda(i) = (i-1)*lambdaint+lambdamin;
    vdata(i,1:2) = v_crit(df,lambda(i),tr,wb);
end
figure(2)
plot(lambda*180/pi,vdata(:,1),'r',lambda*180/pi,vdata(:,2),'b')
xlabel('Head Tube Angle [deg]')
ylabel('Velocity [m/s]')
legend('Weave Critical Velocity','Capsize Critical Velocity')
%---------Vary Trail
clear
df      = 0.6985;                        % front wheel diameter [m]
lambda  = 73*pi/180;                     % head tube angle [rad]
wb      = 1.0287;                        % wheel base [m]
trmin = 0;
trmax = 0.25;
trint = 0.01;
n = (trmax-trmin)/trint+1;
for i = 1:n
    tr(i) = (i-1)*trint+trmin;
    vdata(i,1:2) = v_crit(df,lambda,tr(i),wb);
end
figure(3)
plot(tr,vdata(:,1),'r',tr,vdata(:,2),'b')
xlabel('Trail [m]')
ylabel('Velocity [m/s]')
legend('Weave Critical Velocity','Capsize Critical Velocity')
%---------Vary Wheelbase
clear
df      = 0.6985;                        % front wheel diameter [m]
lambda  = 73*pi/180;                     % head tube angle [rad]
tr      = 0.0602;                        % trail [m]
wbmin = 0.82;
wbmax = 1.36;
wbint = 0.02;
n = ((wbmax-wbmin)/wbint)+1;
for i = 1:n
    wb(i) = (i-1)*wbint+wbmin;
    vdata(i,1:2) = v_crit(df,lambda,tr,wb(i));
end
figure(4)
plot(wb,vdata(:,1),'r',wb,vdata(:,2),'b')
xlabel('Wheelbase [m]')
ylabel('Velocity [m/s]')
legend('Weave Critical Velocity','Capsize Critical Velocity')