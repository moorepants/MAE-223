%    FILE: FINDQ8.AL
%    DATE: MAY 22, 2006
%  AUTHOR: JASON MOORE 
% PROBLEM: PROGRAM TO FIND AN EXPRESSION FOR THE DEPENDENT GENERALIZED
%          COORDINATE (Q8 - FRAME PITCH) OF A FOUR BODY BICYCLE MODEL
%---------------------------------------------------------------------%
%         DEFAULT SETTINGS
%---------------------------------------------------------------------%
OVERWRITE ALL
AUTORHS ALL
%---------------------------------------------------------------------%
%         NEWTONIAN, BODIES, FRAMES, PARTICLES, POINTS
%---------------------------------------------------------------------%
% DECLARE THE INERTIAL REFERENCE FRAME
NEWTONIAN N
% DECLARE THREE INTERMEDIATE FRAMES
% A: YAW FRAME
% B: ROLL FRAME
% E: HEAD TUBE ANGLE FRAME
FRAMES A,B,E
% DECLARE FOUR BODIES
% C: REAR WHEEL
% D: FRAME (INCLUDES RIDER)
% F: FORK
% G: FRONT WHEEL
BODIES C,D,F,G
% DECLARE FOUR POINTS
% NC: REAR CONTACT POINT ON GROUND
% CN: REAR CONTACT POINT ON WHEEL
% NG: FRONT CONTACT POINT ON GROUND
% GN: FRONT CONTACT POINT ON WHEEL
POINTS NC,CN,NG,GN
%---------------------------------------------------------------------%
%         CONSTANTS AND VARIABLES
%---------------------------------------------------------------------%
% RF: RADIUS OF FRONT WHEEL
% RR: RADIUS OF REAR WHEEL
% L1: THE HORIZONTAL DISTANCE FROM THE CENTER OF THE REAR WHEEL TO THE
%     INTERSECTION OF THE HEAD TUBE AXIS
% L2: THE DISTANCE ALONG THE HEAD TUBE AXIS FROM THE PREVIOUS
%     INTERSECTION POINT TO THE INTERSECTION OF THE L3 LINE
% L3: PERPENDICULAR DISTANCE FROM THE HEAD TUBE AXIS TO THE CENTER OF
%     THE FRONT WHEEL
% THETA: COMPLEMENT TO THE HEAD TUBE ANGLE
CONSTANTS RF,RR,L{3},THETA
% DECLARE THE GENERALIZED COORDINATES
% Q1: PERPENDICULAR DISTANCE FROM THE N2> AXIS TO THE REAR CONTACT
%     POINT IN THE GROUND PLANE
% Q2: PERPENDICULAR DISTANCE FROM THE N1> AXIS TO THE REAR CONTACT
%     POINT IN THE GROUND PLANE
% Q3: FRAME YAW ANGLE
% Q4: FRAME ROLL ANGLE
% Q5: REAR WHEEL ROTATION ANGLE
% Q6: STEERING ROTATION ANGLE
% Q7: FRONT WHEEL ROTATION ANGLE
% Q8: FRAME PITCH ANGLE
VARIABLES Q{8}''
%---------------------------------------------------------------------%
%         GENERALIZED SPEEDS
%---------------------------------------------------------------------%
VARIABLES U4',U6'
%---------------------------------------------------------------------%
%         ANGULAR RELATIONSHIPS                                       %
%---------------------------------------------------------------------%
% FRAME YAW
SIMPROT(N,A,3,Q3)
% FRAME ROLL
SIMPROT(A,B,-1,Q4)
% REAR WHEEL ROTATION
SIMPROT(B,C,2,Q5)
% FRAME PITCH
SIMPROT(B,D,-2,Q8)
% HEAD TUBE ANGLE
SIMPROT(D,E,-2,THETA) 
% STEERING ANGLE
SIMPROT(E,F,3,Q6)
% FRONT WHEEL ROTATION
SIMPROT(F,G,2,Q7)
%---------------------------------------------------------------------%
%         POSITION VECTORS
%---------------------------------------------------------------------%
% BUILD THE POSITION VECTOR FROM THE INERTIAL ORIGIN TO NC
P_NO_NC>=Q1*N1>+Q2*N2>
% BUILD THE POSITION VECTOR FROM NC TO CN
P_NC_CN>=0>
% BUILD THE POSITION VECTOR FROM CN TO CO
P_CN_CO>=RR*B3>
% BUILD THE POSITION VECTOR FROM CO TO GO
P_CO_GO>=L1*D1>-L2*F3>+L3*F1>
% BUILD THE POSITION VECTOR FROM GO TO GN
P_GO_GN>=RF*UNITVEC(DOT(F2>,N3>)*F2>-N3>)
% BUILD THE POSITION VECTOR FROM GN TO NG
P_GN_NG>=0>
%---------------------------------------------------------------------%
%        SOLVE FOR Q8, Q8', AND Q8'' IN TERMS OF THE DEPENDENT GENERALIZED COORDINATES
%---------------------------------------------------------------------%
% SET THE N3> COMPONENT OF P_NC_NG> EQUAL TO ZERO
ZERO=DOT(P_NC_NG>,N3>)
% LINEARIZE THE RIGHT HAND SIDE WRT Q8
ZEROL=TAYLOR(ZERO,0:1,Q8=0)
% SOLVE FOR Q8
SOLVE([ZEROL],Q8)
% DEFINE GENERALIZED SPEEDS
Q4'=U4
Q6'=U6
% COMPUTE THE FIRST OF Q8
Q8'=DT(Q8)
%---------------------------------------------------------------------%
%        CHECK EQUATION
%---------------------------------------------------------------------%
CHK1=EVALUATE(Q8,Q6=0)
CHK2=EVALUATE(Q8,Q4=0)
CHK3=EVALUATE(Q8,Q4=0,Q6=0)
%---------------------------------------------------------------------%
%        SAVE OUTPUT
%---------------------------------------------------------------------%
SAVE FINDQ8.ALL