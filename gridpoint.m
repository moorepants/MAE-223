function grid = gridpoint(geometry)
warning off MATLAB:fzero:UndeterminedSyntax
%-------------------Frame Dimensions
alphah  = geometry(1);  % hunch angle [rad]
alphast = geometry(2);  % seat tube angle [rad]
fo      = geometry(3);  % fork offset [m]
hb      = geometry(4);  % bottom bracket height [m]
hs      = geometry(5);  % stem height [m]
lambda  = geometry(6);  % head tube angle [rad]
lc      = geometry(7);  % calf length [m]
lcs     = geometry(8);  % chain stay length [m]
lf      = geometry(9);  % fork length [m]
lhb     = geometry(10); % handle bar length [m]
lla     = geometry(11); % lower arm length [m]
lsp     = geometry(12); % seat post length [m]
lst     = geometry(13); % seat tube length [m]
lth     = geometry(14); % thigh length [m]
lto     = geometry(15); % torso length [m]
lua     = geometry(16); % upper arm length [m]
rf      = geometry(17); % front wheel radius [m]
rh      = geometry(18); % head radius [m]
rr      = geometry(19); % rear wheel radius [m]
wb      = geometry(20); % wheel base [m]
wf      = geometry(21); % front hub width [m]
whb     = geometry(22); % handle bar width [m]
wr      = geometry(23); % rear hub width [m]
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
d = fzero(@leg,[0 0.25],[],x5,x13,z5,z13,lth,lc);
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
d = fzero(@arm,0.25,[],x19,x21,z19,z21,lla,lua);
c=sqrt(lla^2-d^2);
grid(23,1:3) = [grid(19,1)-d,              grid(17,2), grid(19,3)+c              ];
grid(24,1:3) = [grid(23,1),                grid(18,2), grid(23,3)                ];
grid(25,1:3) = [grid(15,1)+rh*cos(alphah), 0,          grid(15,3)+rh*sin(alphah) ];

function zero = leg(d,x5,x13,z5,z13,lth,lc)
zero=lth^2-lc^2-(z13-z5)^2-(x5-x13)^2+2*(z13-z5)*(lc^2-d.^2).^0.5-2*(x5-x13).*d;

function zero = arm(d,x19,x21,z19,z21,lla,lua)
zero=(z21-z19)^2+(x19-x21)^2+lla^2-2*(z21-z19)*(lla^2-d.^2).^0.5-2*(x19-x21)*d-lua^2;