% File: bike_inertia.m
% Date: June 11, 2006
% Author: Jason Moore
% Description:
clear
close all
clc
%-------------------Frame Dimensions
df      = 0.7;                           % front wheel diameter [m]
dr      = 0.6;                           % rear wheel diameter [m]
rf      = df/2;                          % front wheel radius [m]
rr      = dr/2;                          % rear wheel radius [m]
wb      = 1.02;                          % wheel base [m]
hb      = 0.28;                          % bottom bracket height [m]
lcs     = rr+0.1;                        % chain stay length [m]
alphast = 73*pi/180;                     % seat tube angle [rad]
lst     = 0.59;                          % seat tube length [m]
tr      = 0.08;                          % trail [m]
lambda  = atan(3);                       % head tube angle [rad]
fo      = rf*cos(lambda)-tr*sin(lambda); % fork offset [m]
lf      = rf + 0.05;                     % fork length [m]
lsp     = 0.19;                          % seat post length [m]
wr      = 0.127;                         % rear hub width [m]
wf      = 0.098;                         % front hub width [m]
whb     = 0.381;                         % handle bar width [m]
lhb     = 0.1524;                        % handle bar length [m]
hs      = 0.127;                         % stem height [m]
rho     = 7850;                          % frame density [kg/m^3]
%-------------------Human Dimensions
lth     = 0.46;                          % thigh length [m]
lc      = 0.46;                          % calf length [m]
lto     = 0.48;                          % torso length [m]
alphah  = 60*pi/180;                     % hunch angle [rad]
lua     = 0.2794;                        % upper arm length [m]
lla     = 0.3302;                        % lower arm length [m]
rh      = 0.1;                           % head radius [m]
mr      = 75;                            % mass of rider (head, torso, legs, arms) [kg]
%-------------------Grid Point Matrix
% calculates important points from the given frame and human dimensions
grid(1,1:3)  = [0,                                                     0,     0                                                   ];
grid(2,1:3)  = [0,                                                     0,     rr                                                  ];
grid(3,1:3)  = [0,                                                     wr/2,  grid(2,3)                                           ];
grid(4,1:3)  = [0,                                                    -wr/2,  grid(2,3)                                           ];
grid(5,1:3)  = [sqrt(lcs^2-(rr-hb)^2),                                 0,     hb                                                  ];
grid(6,1:3)  = [wb,                                                    0,     0                                                   ];
grid(7,1:3)  = [grid(6,1),                                             0,     rf                                                  ];
grid(8,1:3)  = [grid(6,1),                                             wf/2,  grid(7,3)                                           ];
grid(9,1:3)  = [grid(6,1),                                            -wf/2,  grid(7,3)                                           ];
grid(10,1:3) = [grid(5,1)-lst*cos(alphast),                            0,     grid(5,3)+lst*sin(alphast)                          ];
grid(11,1:3) = [grid(7,1)-fo*sin(lambda)-sqrt(lf^2-fo^2)*cos(lambda),  0,     grid(7,3)-fo*cos(lambda)+sqrt(lf^2-fo^2)*sin(lambda)];
grid(12,1:3) = [grid(11,1)+(grid(11,3)-grid(10,3))/tan(lambda),        0,     grid(10,3)                                          ];
grid(13,1:3) = [grid(10,1)-lsp*cos(alphast),                           0,     grid(10,3)+lsp*sin(alphast)                         ];
% find the zero of the nonlinear relationship for leg position
x5 = grid(5,1); x13 = grid(13,1); z5 = grid(5,3); z13 = grid(13,3);
d=0:0.0001:lc;
zero=lth^2-lc^2-(z13-z5)^2-(x5-x13)^2+2*(z13-z5)*(lc^2-d.^2).^0.5-2*(x5-x13).*d;
for i=1:length(zero)
    if zero(i)<= 0
        index=i;
        break
    end
end
%plot(d,zero,'b',d,zeros(1,length(d)),'r')
d=mean([d(index-1), d(index)]);
b=sqrt(lc^2-d^2);
grid(14,1:3) = [grid(5,1)+d,                 0,           grid(5,3)+b               ];
grid(15,1:3) = [grid(13,1)+lto*cos(alphah),  0,           grid(13,3)+lto*sin(alphah)];
grid(16,1:3) = [grid(12,1)-hs*cos(lambda),   0,           grid(12,3)+hs*sin(lambda) ];
grid(17,1:3) = [grid(16,1),                  whb/2,       grid(16,3)                ];
grid(18,1:3) = [grid(16,1),                 -whb/2,       grid(16,3)                ];
grid(19,1:3) = [grid(17,1)-lhb,              grid(17,2),  grid(17,3)                ];
grid(20,1:3) = [grid(18,1)-lhb,              grid(18,2),  grid(18,3)                ];
grid(21,1:3) = [grid(15,1),                  grid(17,2),  grid(15,3)                ];
grid(22,1:3) = [grid(15,1),                  grid(18,2),  grid(15,3)                ];
% find the zero of the nonlinear relationship for arm position
z21 = grid(21,3); z19 = grid(19,3); x19 = grid(19,1); x21 = grid(21,1);
d=0:0.0001:lla;
zero=(z21-z19)^2+(x19-x21)^2+lla^2-2*(z21-z19)*(lla^2-d.^2).^0.5-2*(x19-x21)*d-lua^2;
% plot(d,zero,'b',d,zeros(1,length(d)),'r')
for i=1:length(zero)
    if zero(i)>= 0
        index=i;
        break
    end
end
d=mean([d(index-1), d(index)]);
c=sqrt(lla^2-d^2);
grid(23,1:3) = [grid(19,1)-d,              grid(17,2), grid(19,3)+c              ];
grid(24,1:3) = [grid(23,1),                grid(18,2), grid(23,3)                ];
grid(25,1:3) = [grid(15,1)+rh*cos(alphah), 0,          grid(15,3)+rh*sin(alphah) ];
%-------------------Element Data Matrix
% column 1: grid point (center for wheels and spheres), starting point for
% others
% column 2: grid point (center for wheels and spheres), ending point for
% others
% column 3: element description
% column 4: element type
% column 5: rigid body element belongs to
% column 6: outer radius for rings and spheres [m], outer radius for tubes [m],
% outer radius for cylinders [m], thickness for rect. prisms [m]
% column 7: tube radius for rings [m], wall thickness for tubes [m], width for
% rect. prisms [m]
ele( 1,1:7) = {2,   2, 'rear wheel',       'ring',     'rwheel', rr,     0.0215 };
ele( 2,1:7) = {3,   5, 'left chain stay',  'tube',     'frame',  0.009,  0.0014 };
ele( 3,1:7) = {4,   5, 'right chain stay', 'tube',     'frame',  0.009,  0.0014 };
ele( 4,1:7) = {3,  10, 'left seat stay',   'tube',     'frame',  0.007,  0.0014 };
ele( 5,1:7) = {4,  10, 'right seat stay',  'tube',     'frame',  0.007,  0.0014 };
ele( 6,1:7) = {5,  13, 'seat tube',        'tube',     'frame',  0.014,  0.0014 };
ele( 7,1:7) = {5,  11, 'down tube',        'tube',     'frame',  0.014,  0.0014 };
ele( 8,1:7) = {10, 12, 'top tube',         'tube',     'frame',  0.0125, 0.0014 };
ele( 9,1:7) = {11, 12, 'head tube',        'tube',     'frame',  0.016,  0.0014 };
ele(10,1:7) = {11, 16, 'fork tube',        'tube',     'fork',   0.0125, 0.0014 };
ele(11,1:7) = {8,  11, 'left fork blade',  'tube',     'fork',   0.011,  0.0014 };
ele(12,1:7) = {9,  11, 'right fork blade', 'tube',     'fork',   0.011,  0.0014 };
ele(13,1:7) = {7,   7, 'front wheel',      'ring',     'fwheel', rf,     0.0215 };
ele(14,1:7) = {17, 18, 'handle bar',       'tube',     'fork',   0.011,  0.0014 };
ele(15,1:7) = {17, 19, 'left handle',      'tube',     'fork',   0.011,  0.0014 };
ele(16,1:7) = {18, 20, 'right handle',     'tube',     'fork',   0.011,  0.0014 };
ele(17,1:7) = {5,  14, 'calfs',            'rprism',   'frame',  0.13,   0.40   };
ele(18,1:7) = {13, 14, 'thighs',           'rprism',   'frame',  0.15,   0.40   };
ele(19,1:7) = {13, 15, 'torso',            'rprism',   'frame',  0.23,   0.40   };
ele(20,1:7) = {25, 25, 'head',             'sphere',   'frame',  rh,     0      };
ele(21,1:7) = {21, 23, 'left upper arm',   'cylinder', 'frame',  0.1,    0      };
ele(22,1:7) = {22, 24, 'right upper arm',  'cylinder', 'frame',  0.1,    0      };
ele(23,1:7) = {19, 23, 'left lower arm',   'cylinder', 'frame',  0.1,    0      };
ele(24,1:7) = {20, 24, 'right lower arm',  'cylinder', 'frame',  0.1,    0      };
%-------------------Draw 2D Bicycle
figure(1)
axis([-rr-0.2 wb+rf+0.8 0 wb+rf+0.8])
hold on
for i = 1:size(ele,1)  
    switch ele{i,4}
        case {'ring','sphere'} % plot circles
            center = [grid(ele{i,1},1);grid(ele{i,1},3)];
            radius = ele{i,6};
            thetaplot = 0:0.001:2*pi;
            xplot = radius*cos(thetaplot)+center(1);
            yplot = radius*sin(thetaplot)+center(2);            
            plot(xplot,yplot)
        case {'tube','rprism','cylinder'} % plot lines
            if grid(ele{i,1},1) < grid(ele{i,2},1)
                x1 = grid(ele{i,1},1);
                x2 = grid(ele{i,2},1);
                z1 = grid(ele{i,1},3);
                z2 = grid(ele{i,2},3);
                xplot=x1:0.001:x2;
                yplot=(z2-z1)/(x2-x1).*(xplot-x1)+z1;
                plot(xplot,yplot)
            else
                x2 = grid(ele{i,1},1);
                x1 = grid(ele{i,2},1);
                z2 = grid(ele{i,1},3);
                z1 = grid(ele{i,2},3);
                xplot=x1:0.001:x2;
                yplot=(z2-z1)/(x2-x1).*(xplot-x1)+z1;
                plot(xplot,yplot)                
            end                
    end
end
hold off
%-------------------Calculate Element Length
for i = 1:size(ele,1)
    x1 = grid(ele{i,1},1);
    x2 = grid(ele{i,2},1);
    y1 = grid(ele{i,1},2);
    y2 = grid(ele{i,2},2);
    z1 = grid(ele{i,1},3);
    z2 = grid(ele{i,2},3);
    elelength(i)=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
end
%-------------------Calculate Element Mass
elemass = zeros(size(ele,1),1); % build empty vector
% define masses of human body parts and wheels
elemass(1)  = 6.67*rr;
elemass(13) = 8.57*rf;
elemass(17) = .15*mr;
elemass(18) = .15*mr;
elemass(19) = .3*mr;
elemass(20) = .1*mr;
elemass(21) = .075*mr;
elemass(22) = .075*mr;
elemass(23) = .075*mr;
elemass(24) = .075*mr;
% calculate the mass of the frame and fork tube elements
for i = 1:size(ele,1)
    switch ele{i,4}
        case 'tube'
        elemass(i) = rho*pi*elelength(i)*(2*ele{i,6}*ele{i,7}-(ele{i,7})^2);
    end
end
%-------------------Calculate Total Mass of Each Rigid Body
% initialize mass summation variables
framemass = 0;
forkmass = 0;
rwheelmass = 0;
fwheelmass = 0;
% sum the mass from each element for each body
for i = 1:size(ele,1)
    switch ele{i,5}
        case 'frame'
            framemass = framemass + elemass(i);
        case 'fork'
            forkmass = forkmass + elemass(i);
        case 'rwheel'
            rwheelmass = rwheelmass + elemass(i);
        case 'fwheel'
            fwheelmass = fwheelmass + elemass(i);
    end
end
%-------------------Element Center of Mass
hold on
for i = 1:size(ele,1)
    x1 = grid(ele{i,1},1);
    x2 = grid(ele{i,2},1);
    y1 = grid(ele{i,1},2);
    y2 = grid(ele{i,2},2);
    z1 = grid(ele{i,1},3);
    z2 = grid(ele{i,2},3);
    elecenter(i,1:3) = [(x1+x2)/2,(y1+y2)/2,(z1+z2)/2]; % midpoint formula
    plot(elecenter(i,1),elecenter(i,3),'ro') % plot element centers on 2D figure
end
hold off
%-------------------Rigid Body Center of Mass
% initialize center of mass vectors
framecenter = [0,0,0];
forkcenter = [0,0,0];
rwheelcenter = [0,0,0];
fwheelcenter = [0,0,0];
for i = 1:size(ele,1)
    switch ele{i,5}
        case 'frame'
            framecenter = framecenter + [elecenter(i,1).*elemass(i),elecenter(i,2).*elemass(i),elecenter(i,3).*elemass(i)];
        case 'fork'
            forkcenter = forkcenter + [elecenter(i,1).*elemass(i),elecenter(i,2).*elemass(i),elecenter(i,3).*elemass(i)];
        case 'rwheel'
            rwheelcenter = rwheelcenter + [elecenter(i,1).*elemass(i),elecenter(i,2).*elemass(i),elecenter(i,3).*elemass(i)];
        case 'fwheel'
            fwheelcenter = fwheelcenter + [elecenter(i,1).*elemass(i),elecenter(i,2).*elemass(i),elecenter(i,3).*elemass(i)];
    end
end
% divide sum by total mass of each body
framecenter = framecenter./framemass;
forkcenter = forkcenter./forkmass;
rwheelcenter = rwheelcenter./rwheelmass;
fwheelcenter = fwheelcenter./fwheelmass;
hold on
plot([framecenter(1);forkcenter(1);rwheelcenter(1);fwheelcenter(1)],[framecenter(3);forkcenter(3);rwheelcenter(3);fwheelcenter(3)],'go')
hold off
%-------------------Element Inertia Tensor
% about local reference frame and element center of mass
for i = 1:size(ele,1)
    switch ele{i,4}
        case 'tube' % z axis is always along length of tube
            eleinertialocal{i} = [1/12*elemass(i)*(3*((ele{i,6})^2+(ele{i,6}-ele{i,7})^2)+(elelength(i))^2),0,0;
                                  0,1/12*elemass(i)*(3*((ele{i,6})^2+(ele{i,6}-ele{i,7})^2)+(elelength(i))^2),0;
                                  0,0,1/2*elemass(i)*((ele{i,6})^2+(ele{i,6}-ele{i,7})^2)];
        case 'rprism' % z axis is always along length of prism
            eleinertialocal{i} = [1/12*elemass(i)*((ele{i,7})^2+(elelength(i))^2),0,0;
                                  0,1/12*elemass(i)*((ele{i,6})^2+(elelength(i))^2),0;
                                  0,0,1/12*elemass(i)*((ele{i,6})^2+(ele{i,7})^2)];
        case 'cylinder'
            eleinertialocal{i} = [1/12*elemass(i)*(3*(ele{i,6})^2+(elelength(i))^2),0,0;
                                  0,1/12*elemass(i)*(3*(ele{i,6})^2+(elelength(i))^2),0;
                                  0,0,1/2*elemass(i)*(ele{i,6})^2];
        case 'ring'
            eleinertialocal{i} = [1/8*(5*(ele{i,7})^2+4*(ele{i,6}-ele{i,7})^2)*elemass(i),0,0;
                                  0, (3/4*(ele{i,7})^2+(ele{i,6}-ele{i,7})^2)*elemass(i),0;
                                  0,0,1/8*(5*(ele{i,7})^2+4*(ele{i,6}-ele{i,7})^2)*elemass(i)];
        case 'sphere'
            eleinertialocal{i} = [2/5*elemass(i)*(ele{i,6})^2,0,0;
                                  0,2/5*elemass(i)*(ele{i,6})^2,0;
                                  0,0,2/5*elemass(i)*(ele{i,6})^2];
    end
end
%-------------------Element Unit Vectors
%local z vector
for i = 1:size(ele,1)
    switch ele{i,4}
        case {'ring','sphere'}
            uveczlocal(i,1:3) = [0,0,1];
        otherwise
            veczlocal(i,1:3) = [grid(ele{i,2},1)-grid(ele{i,1},1),grid(ele{i,2},2)-grid(ele{i,1},2),grid(ele{i,2},3)-grid(ele{i,1},3)];
            mag(i,1) = sqrt(sum((veczlocal(i,1:3).^2)'));
            uveczlocal(i,1:3) = veczlocal(i,1:3)./mag(i);
    end
end
%local y vector
for i = 1:size(ele,1)
    switch ele{i,3}
        case {'right seat stay','right chain stay'}
            vecylocal(i,1:3) = cross(uveczlocal(5,1:3),uveczlocal(3,1:3));
            mag(i,1) = sqrt(sum((vecylocal(i,1:3).^2)'));
            uvecylocal(i,1:3) = vecylocal(i,1:3)./mag(i);
        case {'left seat stay','left chain stay'}
            vecylocal(i,1:3) = cross(uveczlocal(4,1:3),uveczlocal(2,1:3));
            mag(i,1) = sqrt(sum((vecylocal(i,1:3).^2)'));
            uvecylocal(i,1:3) = vecylocal(i,1:3)./mag(i);
        case {'right fork blade','left fork blade'}
            vecylocal(i,1:3) = cross(uveczlocal(11,1:3),uveczlocal(12,1:3));
            mag(i,1) = sqrt(sum((vecylocal(i,1:3).^2)'));
            uvecylocal(i,1:3) = vecylocal(i,1:3)./mag(i);
        case 'handle bar'
            uvecylocal(i,1:3) = [1,0,0];
        otherwise
            uvecylocal(i,1:3) = [0,1,0];
    end
end
%local x vector
for i = 1:size(ele,1)
    uvecxlocal(i,1:3) = cross(uvecylocal(i,1:3),uveczlocal(i,1:3));
end
%-------------------Element Direction Cosines Relative to Global Frame
xhat = [1,0,0];
yhat = [0,1,0];
zhat = [0,0,1];
for i = 1:size(ele,1)
    dircos{i} = [dot(xhat,uvecxlocal(i,1:3)),dot(xhat,uvecylocal(i,1:3)),dot(xhat,uveczlocal(i,1:3));
                 dot(yhat,uvecxlocal(i,1:3)),dot(yhat,uvecylocal(i,1:3)),dot(yhat,uveczlocal(i,1:3));
                 dot(zhat,uvecxlocal(i,1:3)),dot(zhat,uvecylocal(i,1:3)),dot(zhat,uveczlocal(i,1:3))];
end
%-------------------Rotate Element Interia Tensors to Global Frame
for i = 1:size(ele,1)
    eleinertiaglobal{i} = dircos{i}*eleinertialocal{i}*(dircos{i})';
end
%-------------------Translate Inertia Tensors to the Center of Mass of Each Body
for i = 1:size(ele,1)
    switch ele{i,5}
        case 'frame'
            a = framecenter(1) - elecenter(i,1);
            b = framecenter(2) - elecenter(i,2);
            c = framecenter(3) - elecenter(i,3);
            eleinertiatrans{i} = eleinertiaglobal{i} + elemass(i)*[b^2+c^2, -a*b, -a*c;-a*b,c^2+a^2,-b*c;-a*c,-b*c,a^2+b^2];
        case 'fork'
            a = forkcenter(1) - elecenter(i,1);
            b = forkcenter(2) - elecenter(i,2);
            c = forkcenter(3) - elecenter(i,3);
            eleinertiatrans{i} = eleinertiaglobal{i} + elemass(i)*[b^2+c^2, -a*b, -a*c;-a*b,c^2+a^2,-b*c;-a*c,-b*c,a^2+b^2];
        case 'rwheel'
            a = rwheelcenter(1) - elecenter(i,1);
            b = rwheelcenter(2) - elecenter(i,2);
            c = rwheelcenter(3) - elecenter(i,3);
            eleinertiatrans{i} = eleinertiaglobal{i} + elemass(i)*[b^2+c^2, -a*b, -a*c;-a*b,c^2+a^2,-b*c;-a*c,-b*c,a^2+b^2];
        case 'fwheel'
            a = fwheelcenter(1) - elecenter(i,1);
            b = fwheelcenter(2) - elecenter(i,2);
            c = fwheelcenter(3) - elecenter(i,3);
            eleinertiatrans{i} = eleinertiaglobal{i} + elemass(i)*[b^2+c^2, -a*b, -a*c;-a*b,c^2+a^2,-b*c;-a*c,-b*c,a^2+b^2];
    end
end
%sum inertia tensors for each rigid body
frameinertia = zeros(3);
forkinertia = zeros(3);
rwheelinertia = zeros(3);
fwheelinertia = zeros(3);
for i = 1:size(ele,1)
    switch ele{i,5}
        case 'frame'
            frameinertia = frameinertia + eleinertiatrans{i};
        case 'fork'
            forkinertia = forkinertia + eleinertiatrans{i};
        case 'rwheel'
            rwheelinertia = rwheelinertia + eleinertiatrans{i};
        case 'fwheel'
            fwheelinertia = fwheelinertia + eleinertiatrans{i};
    end
end
% rotate fork inertia tensor into fork reference frame
theta = pi/2 - lambda;
Cfork = [ cos(theta), 0, sin(theta);
          0,          1,          0;
         -sin(theta), 0, cos(theta) ];
forkinertia = Cfork*forkinertia*Cfork';
% Autolev constants
RF = rf;
RR = rr;
THETA = theta;
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
MC = rwheelmass;
MD = framemass;
MF = forkmass;
MG = fwheelmass;
IC11 = rwheelinertia(1,1);
IC22 = rwheelinertia(2,2);
IC33 = rwheelinertia(3,3);
IC12 = rwheelinertia(1,2);
IC23 = rwheelinertia(2,3);
IC31 = rwheelinertia(3,1);
ID11 = frameinertia(1,1);
ID22 = frameinertia(2,2);
ID33 = frameinertia(3,3);
ID12 = frameinertia(1,2);
ID23 = frameinertia(2,3);
ID31 = frameinertia(3,1);
IF11 = forkinertia(1,1);
IF22 = forkinertia(2,2);
IF33 = forkinertia(3,3);
IF12 = forkinertia(1,2);
IF23 = forkinertia(2,3);
IF31 = forkinertia(3,1);
IG11 = fwheelinertia(1,1);
IG22 = fwheelinertia(2,2);
IG33 = fwheelinertia(3,3);
IG12 = fwheelinertia(1,2);
IG23 = fwheelinertia(2,3);
IG31 = fwheelinertia(3,1);