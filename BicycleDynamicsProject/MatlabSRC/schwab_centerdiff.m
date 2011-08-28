% File: schwab.m
% Date: March 14, 2007
% Author: Jason Moore
% Description: Compares the bicycle model equations of motion from
% Schwab to the ones developed by Moore in Autolev. Updated to
% difference method and different delta.
clear
close all
clc
%---------Define Schwab's Parameters
%----General
w     = 1.02;            % wheelbase [m]
t     = 0.08;            % trail [m]
alpha = atan(3);         % head tube angle [m]
g     = 9.81;            % acceleration due to gravity [m/s^2]
%----Rear Wheel
Rrw   = 0.3;             % radius [m]
xrw   = 0;               % mass center x coordinate [m]
yrw   = 0;               % mass center y coordinate [m]
zrw   = -Rrw;            % mass center z coordinate [m]
cmrw  = [xrw;yrw;zrw];   % mass center vector [m]
mrw   = 2;               % mass [kg]
Axx   = 0.06;            % moment of inertia about x axis [kg*m^2]
Ayy   = 0.12;            % moment of inertia about y axis [kg*m^2]
Azz   = 0.06;            % moment of inertia about z axis [kg*m^2]
Ai    = [Axx,   0,   0;
           0, Ayy,   0;
           0,   0, Azz]; % inertia tensor
%----Rear Frame
xrf   = 0.3;             % mass center x coordinate [m]
yrf   = 0;               % mass center y coordinate [m]
zrf   = -0.9;            % mass center z coordinate [m]
cmrf  = [xrf;yrf;zrf];   % mass center vector [m]
mrf   = 85;              % mass [kg]
Bxx   = 9.2;             % moment of inertia about x axis [kg*m^2]
Bxy   = 0;               % product of inertia [kg*m^2]
Bxz   = 2.4;             % product of inertia [kg*m^2]
Byy   = 11;              % moment of inertia about y axis [kg*m^2]
Byz   = 0;               % product of inertia [kg*m^2]
Bzz   = 2.8;             % moment of inertia about z axis [kg*m^2]
Bi    = [Bxx, Bxy, Bxz;
         Bxy, Byy, Byz;
         Bxz, Byz, Bzz]; % inertia tensor
%----Front Frame
xff   = 0.9;             % mass center x coordinate [m]
yff   = 0;               % mass center y coordinate [m]
zff   = -0.7;            % mass center z coordinate [m]
cmff  = [xff;yff;zff];   % mass center vector [m]
mff   = 4;               % mass [kg]
Cxx   = 0.0546;          % moment of inertia about x axis [kg*m^2]
Cxy   = 0;               % product of inertia [kg*m^2]
Cxz   = -0.0162;         % product of inertia [kg*m^2]
Cyy   = 0.06;            % moment of inertia about y axis [kg*m^2]
Cyz   = 0;               % product of inertia [kg*m^2]
Czz   = 0.0114;          % moment of inertia about z axis [kg*m^2]
Ci    = [Cxx, Cxy, Cxz;
         Cxy, Cyy, Cyz;
         Cxz, Cyz, Czz]; % inertia tensor
%----Front Wheel
Rfw   = 0.35;            % radius [m]
xfw   = w;               % mass center x coordinate [m]
yfw   = 0;               % mass center y coordinate [m]
zfw   = -Rfw;            % mass center z coordinate [m]
cmfw  = [xfw;yfw;zfw];   % mass center vector [m]
mfw   = 3;               % mass [kg]
Dxx   = 0.14;            % moment of inertia about x axis [kg*m^2]
Dyy   = 0.28;            % moment of inertia about y axis [kg*m^2]
Dzz   = 0.14;            % moment of inertia about z axis [kg*m^2]
Di    = [Dxx,   0,   0;
           0, Dyy,   0;
           0,   0, Dzz]; % inertia tensor
%---------Convert Schwab's Parameters To Moore's Parameters
wb = w;         % wheelbase [m]
tr = t;         % trail [m]
lambda = alpha; % head tube angle [m]
rr = Rrw;       % rear wheel radius [m]
rf = Rfw;       % front wheel radius [m]
dr = rr*2;      % rear wheel diamter [m]
df = rf*2;      % front wheel diameter [m]
fo = rf*cos(lambda)-tr*sin(lambda); % fork offset [m]
%----rotate to moore's global reference frame
sCm = [1,  0,  0;
       0, -1,  0;
       0,  0, -1]; % direction cosine matrix (Moore relative to Schwab)
rwheelcenter = sCm*cmrw;
framecenter = sCm*cmrf;
forkcenter = sCm*cmff;
fwheelcenter = sCm*cmfw;
rwheelinertia = sCm*Ai*sCm';
frameinertia = sCm*Bi*sCm';
forkinertia = sCm*Ci*sCm';
fwheelinertia = sCm*Di*sCm';
%----rotate fork reference frame through head tube angle
theta = pi/2 - lambda; % complement of the head tube angle [rad]
Cfork = [ cos(theta), 0, sin(theta);
          0,          1,          0;
         -sin(theta), 0, cos(theta) ];
forkinertia = Cfork*forkinertia*Cfork';
%---------Define Autolev Constants
G = g; 
IC11 = rwheelinertia(1,1);
IC22 = rwheelinertia(2,2);
IC33 = rwheelinertia(3,3);
ID11 = frameinertia(1,1);
ID12 = frameinertia(1,2);
ID22 = frameinertia(2,2);
ID23 = frameinertia(2,3);
ID31 = frameinertia(3,1);
ID33 = frameinertia(3,3);
IF11 = forkinertia(1,1);
IF12 = forkinertia(1,2);
IF22 = forkinertia(2,2);
IF23 = forkinertia(2,3);
IF31 = forkinertia(3,1);
IF33 = forkinertia(3,3);
IG11 = fwheelinertia(1,1);
IG22 = fwheelinertia(2,2);
IG33 = fwheelinertia(3,3);
L1 = wb-fo*sin(lambda)-((dr-df)/2+fo*cos(lambda))/tan(lambda);
L2 = ((dr-df)/2+fo*cos(lambda))/sin(lambda);
L3 = fo;
L4 = framecenter(1) - rwheelcenter(1);
L5 = framecenter(3) - rwheelcenter(3);
temp1 = forkcenter(1) - fwheelcenter(1);
temp2 = forkcenter(3) - fwheelcenter(3);
temp3 = temp1^2 + temp2^2;
temp4 = temp1/temp2;
L6 = sqrt((temp3)/(1+(tan(pi/2-theta-atan(temp4)))^2));
L7 = sqrt(temp3-L6^2);
MC = mrw;
MD = mrf;
MF = mff;
MG = mfw;
RF = rf;
RR = rr;
THETA = theta;
autolev_constants = [G;
                     IC11;IC22;IC33;
                     ID11;ID12;ID22;ID23;ID31;ID33;
                     IF11;IF12;IF22;IF23;IF31;IF33;
                     IG11;IG22;IG33;
                     L1;L2;L3;L4;L5;L6;L7;
                     MC;MD;MF;MG;
                     RF;RR;THETA];
%-------------------Construct the Stability Matrix and Calculate
%-------------------Eigenvalues for Various Velocities
delta = 1e-11; % perturbance value
vmax = 10; % max foward velocity of the rear wheel to be calculated
n = 1000; % number of iterations
for i=1:n  
    v(i) = (i-1)/n*vmax; % ith velocity
    NU5(i) = v(i)/RR; % ith angular velocity of the rear wheel 
    % "evalprimes" computes the derivatives of the state variables, the
    % equations of motion were generated in Autolev
    % compute nominal solution for the ith velocity
    nominal(:,i) = evalprimes([0;0;0;0;0;0;0;0;NU5(i);0],autolev_constants);
    % build the stability matrix by numerically calculating the partial
    % derivatives of each differential equation with respect to each state
    % variable
    for j=1:10;
        perturb1 = [0;0;0;0;0;0;0;0;NU5(i);0]; %initialize function input
        perturb2 = [0;0;0;0;0;0;0;0;NU5(i);0]; %initialize function input
        perturb1(j) = perturb1(j) + delta; %perturb the jth variable
        perturb2(j) = perturb2(j) - delta;
        %solve differential equations for perturbance
        prime1 = evalprimes(perturb1,autolev_constants);
        prime2 = evalprimes(perturb2,autolev_constants);
        %compute partial derivative
        matrix(:,j) = (prime1-prime2)./2./delta;
    end
    % reduce stability matrix to 4 x 4 matrix for steer and roll
    stab(1,1:4) = [matrix(4,4)  matrix(4,6)  matrix(4,8)  matrix(4,10) ];
    stab(2,1:4) = [matrix(6,4)  matrix(6,6)  matrix(6,8)  matrix(6,10) ];
    stab(3,1:4) = [matrix(8,4)  matrix(8,6)  matrix(8,8)  matrix(8,10) ];
    stab(4,1:4) = [matrix(10,4) matrix(10,6) matrix(10,8) matrix(10,10)];
    if i == n/2
        stab_moore = stab;
        v_moore = v(i);
    end
    % calculate the eigenvalues for the reduced stability matrix
    eigval(1:4,i)=eig(stab);
end
%-------------------Organize Eigenvalue Matrix
% set first column in organized matrix to the first column of the original
% matrix
eigorg(:,1)=eigval(:,1);
% rearrange the eigenvalue columns by calcuating the absolute value of the
% difference between values of the successicve column and the preceeding
% column
for i = 1:length(v)-1
    for j = 1:4
        first = eigorg(j,i);
        for k =1:4
            second = eigval(k,i+1);
            diff(k) = abs(second - first);
        end
        [mindiff indice]=min(diff);
        eigorg(j,i+1)=eigval(indice,i+1);
    end
end
%-------------------Extract Weave, Capsize, and Caster Eigenmodes
found = 0;
for i = 1:size(eigorg,1)
    if isreal(eigorg(i,:)) == 0 & found==0
        weave(1,:)=real(eigorg(i,:));
        found = 1;
    elseif isreal(eigorg(i,:)) == 0 & found==1
        weave(2,:)=real(eigorg(i,:));
    elseif isreal(eigorg(i,:)) == 1 & ...
           (eigorg(i,size(eigorg,2))-eigorg(i,1))>0
        capsize=real(eigorg(i,:));
    else
        caster=real(eigorg(i,:));
    end
end
%-------------------Plot Eigenmodes vs Velocity
figure(1)
hold on
axis([0 10 -10 10])
title('Eigenmodes')
xlabel('Velocity [m/s]')
ylabel('Eigenvalues (Real) [1/s]')
plot(v,weave(1,:),'.b')
plot(v,capsize,'.y')
plot(v,caster,'.r')
legend('Weave','Capsize','Caster')
plot(v,weave(2,:),'.b')
plot(v,zeros(length(v),1),'k') %plot horizontal line at zero
hold off
%-------------------Calculate Critical Velocities
% find the velocity at which the weave mode becomes stable
for i=1:length(weave(1,:))
    if weave(1,i)<=0
        index=i;
        break
    end
end
vw=mean([v(index-1),v(index)]) % weave critical velocity
% find the velocity at which the capsize mode becomes unstable
for i=1:length(capsize)
    if capsize(i)>=0
        index=i;
        break
    end
end
vc=mean([v(index-1),v(index)]) % capsize critical velocity
%-------------------Plot the Critical Velocities
hold on
plot(vw,0,'ok',zeros(1,20)+vw,-19:1:0,'k',vc,0,'ok',zeros(1,20)+vc,...
    -19:1:0,'k')
hold off
mv = v;
mweave = weave;
mcapsize = capsize;
mcaster = caster;
mvw = vw;
mvc = vc;
%-------------------Convert Schwab's Canonical Form
M  = [  80.81210000000002,   2.32343142623549;
         2.32343142623549,   0.30126570934256];
K0 = [-794.119500000000, -25.739089291258;
       -25.739089291258,  -8.139414705882];
K2 = [0,76.40620875965657;0,2.67560553633218];
C1 = [0,33.77386947593010;-0.84823447825693,1.70696539792387];
eigval = zeros(4,1000);
n=1000;
vmax = 10;
for i=1:n  
    v(i) = (i-1)/n*vmax;
    C = C1.*v(i);
    K = K0 + K2.*v(i)^2;
    stab = zeros(4);
    stab(1,3)=1;
    stab(2,4)=1;
    stab(3,1)=(K(2,1)/M(2,2)-K(1,1)/M(1,2))/(M(1,1)/M(1,2)-M(2,1)/M(2,2));
    stab(3,2)=(K(2,2)/M(2,2)-K(1,2)/M(1,2))/(M(1,1)/M(1,2)-M(2,1)/M(2,2));
    stab(3,3)=(C(2,1)/M(2,2)-C(1,1)/M(1,2))/(M(1,1)/M(1,2)-M(2,1)/M(2,2));
    stab(3,4)=(C(2,2)/M(2,2)-C(1,2)/M(1,2))/(M(1,1)/M(1,2)-M(2,1)/M(2,2));
    stab(4,1)=(K(2,1)/M(2,1)-K(1,1)/M(1,1))/(M(1,2)/M(1,1)-M(2,2)/M(2,1));
    stab(4,2)=(K(2,2)/M(2,1)-K(1,2)/M(1,1))/(M(1,2)/M(1,1)-M(2,2)/M(2,1));
    stab(4,3)=(C(2,1)/M(2,1)-C(1,1)/M(1,1))/(M(1,2)/M(1,1)-M(2,2)/M(2,1));
    stab(4,4)=(C(2,2)/M(2,1)-C(1,2)/M(1,1))/(M(1,2)/M(1,1)-M(2,2)/M(2,1));
    if i == n/2
        stab_schwab = stab;
        v_schwab = v(i);
    end
    eigval(1:4,i)=eig(stab);
end
%-------------------Organize Eigenvalue Matrix
% set first column in organized matrix to the first column of the original
% matrix
eigorg(:,1)=eigval(:,1);
% rearrange the eigenvalue columns by calcuating the absolute value of the
% difference between values of the successicve column and the preceeding
% column
for i = 1:length(v)-1
    for j = 1:4
        first = eigorg(j,i);
        for k =1:4
            second = eigval(k,i+1);
            diff(k) = abs(second - first);
        end
        [mindiff indice]=min(diff);
        eigorg(j,i+1)=eigval(indice,i+1);
    end
end
%-------------------Extract Weave, Capsize, and Caster Eigenmodes
found = 0;
for i = 1:size(eigorg,1)
    if isreal(eigorg(i,:)) == 0 & found==0
        weave(1,:)=real(eigorg(i,:));
        found = 1;
    elseif isreal(eigorg(i,:)) == 0 & found==1
        weave(2,:)=real(eigorg(i,:));
    elseif isreal(eigorg(i,:)) == 1 &...
           (eigorg(i,size(eigorg,2))-eigorg(i,1))>0
        capsize=real(eigorg(i,:));
    else
        caster=real(eigorg(i,:));
    end
end
%-------------------Plot Eigenmodes vs Velocity
figure(2)
hold on
axis([0 10 -10 10])
title('Eigenmodes')
xlabel('Velocity [m/s]')
ylabel('Eigenvalues (Real) [1/s]')
plot(v,weave(1,:),'.b')
plot(v,capsize,'.y')
plot(v,caster,'.r')
legend('Weave','Capsize','Caster')
plot(v,weave(2,:),'.b')
plot(v,zeros(length(v),1),'k') %plot horizontal line at zero
hold off
%-------------------Calculate Critical Velocities
% find the velocity at which the weave mode becomes stable
for i=1:length(weave(1,:))
    if weave(1,i)<=0
        index=i;
        break
    end
end
vw=mean([v(index-1),v(index)]) % weave critical velocity
% find the velocity at which the capsize mode becomes unstable
for i=1:length(capsize)
    if capsize(i)>=0
        index=i;
        break
    end
end
vc=mean([v(index-1),v(index)]) % capsize critical velocity
%-------------------Plot the Critical Velocities
hold on
plot(vw,0,'ok',zeros(1,20)+vw,-19:1:0,'k',vc,0,'ok',zeros(1,20)+vc,...
    -19:1:0,'k')
hold off
sv = v;
sweave = weave;
scapsize = capsize;
scaster = caster;
svw = vw;
svc = vc;
figure(3)
hold on
axis([0 10 -10 10])
title('Eigenmodes')
xlabel('Velocity [m/s]')
ylabel('Eigenvalues (Real) [1/s]')
plot(mv,mweave(1,:),'b',sv,sweave(1,:),'r',mv,mweave(2,:),'b',sv,...
    sweave(2,:),'r')
legend('Autolevs Equations','Schwabs Equations')
plot(mv,mcapsize,'b',mv,mcaster,'b',sv,scapsize,'r',sv,scaster,'r')
plot(sv,zeros(length(v),1),'k') %plot horizontal line at zero
hold off
%-------------------Plot the Critical Velocities
hold on
plot(mvw,0,'ob',zeros(1,20)+mvw,-19:1:0,'b',mvc,0,'ob',zeros(1,20)+mvc,...
    -19:1:0,'b')
plot(svw,0,'or',zeros(1,20)+svw,-19:1:0,'r',svc,0,'or',zeros(1,20)+svc,...
    -19:1:0,'r')
hold off

