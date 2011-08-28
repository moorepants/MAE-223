%---------------------------------------------------------------------%
%    FILE: BICYCLE13.AL
%    DATE: JULY 6, 2006
%  AUTHOR: JASON MOORE
% PROBLEM: GENERATES THE EQUATIONS OF MOTION FOR A FIXED-RIDER 
%          NO-HANDS BICYCLE MADE UP OF FOUR RIGID BODIES.
%---------------------------------------------------------------------%
%         DEFAULT SETTINGS
%---------------------------------------------------------------------%
AUTOZ ON
AUTORHS ON
OVERWRITE ALL
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
% ALL DESCRIPTIONS ARE WITH REFERENCE TO THE NOMINAL UPRIGHT BICYCLE
% POSITION
% RF: RADIUS OF FRONT WHEEL
% RR: RADIUS OF REAR WHEEL
% L1: THE HORIZONTAL DISTANCE FROM THE CENTER OF THE REAR WHEEL TO THE
%     INTERSECTION OF THE HEAD TUBE AXIS
% L2: THE DISTANCE ALONG THE HEAD TUBE AXIS FROM THE PREVIOUS
%     INTERSECTION POINT TO THE INTERSECTION OF THE L3 LINE
% L3: PERPENDICULAR DISTANCE FROM THE HEAD TUBE AXIS TO THE CENTER OF
%     THE FRONT WHEEL
% L4: THE HORIZONTAL DISTANCE FROM THE CENTER OF THE REAR WHEEL TO THE
%     CENTER OF MASS OF THE FRAME/RIDER
% L5: THE VERTICAL DISTANCE FROM THE CENTER OF THE REAR WHEEL TO THE
%     TO THE CENTER OF MASS OF THE FRAME/RIDER
% L6: THE DISTANCE FROM THE CENTER OF THE FRONT WHEEL TO THE LOCATION
%     OF THE CENTER OF MASSS OF THE FRONT FORK THAT IS PERDENDICULAR
%     TO THE HEAD TUBE AXIS
% L7: THE DISTANCE FROM THE CENTER OF THE FRONT WHEEL TO THE LOCATION
%     OF THE CENTER OF MASS OF THE FRONT FORK THAT IS PARALLEL TO THE
%     HEAD TUBE AXIS
% THETA: COMPLEMENT TO THE HEAD TUBE ANGLE
CONSTANTS RF,RR,L{7},THETA
% ACCELERATION DUE TO GRAVITY
CONSTANTS G
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
VARIABLES Q{7}',Q8
%---------------------------------------------------------------------%
%         GENERALIZED SPEEDS
%---------------------------------------------------------------------%
MOTIONVARIABLES' U{7}'
%---------------------------------------------------------------------%
%         MASS AND INERTIA PROPERTIES
%---------------------------------------------------------------------%
MASS C=MC,D=MD,F=MF,G=MG
INERTIA C,IC11,IC22,IC33
INERTIA D,ID11,ID22,ID33,ID12,ID23,ID31
INERTIA F,IF11,IF22,IF33,IF12,IF23,IF31
INERTIA G,IG11,IG22,IG33
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
% PITCH WAS CHOSEN AS THE DEPENDENT GENERALIZED COORDINATE AND SOLVED
% FOR IN TERMS OF ROLL
% (Q4) AND STEER (Q6) IN FINDQ8.AL
Q8=(L2*COS(THETA)*COS(Q4)+L3*(SIN(Q4)*SIN(Q6)-SIN(THETA)*COS(Q4)*COS(&
   Q6))-RR*COS(Q4)-RF*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4)&
   )^2)/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^0.5)/(COS(&
   Q4)*(L1+L2*SIN(THETA)+L3*COS(THETA)*COS(Q6)+2*RF*COS(THETA)*SIN(Q6&
   )*(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))/(1-(SIN(Q4)*COS(Q6)&
   +SIN(THETA)*SIN(Q6)*COS(Q4))^2)^0.5+RF*COS(THETA)*SIN(Q6)*(SIN(Q4)&
   *COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))*(-1+(SIN(Q4)*COS(Q6)+SIN(THET&
   A)*SIN(Q6)*COS(Q4))^2)/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(&
   Q4))^2)^1.5))
SIMPROT(B,D,-2,Q8)
% HEAD TUBE ANGLE
SIMPROT(D,E,-2,THETA)
% STEERING ANGLE
SIMPROT(E,F,3,Q6)
% FRONT WHEEL ROTATION
SIMPROT(F,G,2,Q7)
%---------------------------------------------------------------------%
%         ANGULAR VELOCITIES
%---------------------------------------------------------------------%
W_A_N>=U3*A3>
W_B_A>=-U4*B1>
W_C_B>=U5*C2>
% THE DERIVATIVE OF Q8 WAS CALCULATED IN FINDQ8.AL
W_D_B>=-((L2*COS(THETA)*COS(Q4)+L3*(SIN(Q4)*SIN(Q6)-SIN(THETA)*COS(Q4&
        )*COS(Q6))-RR*COS(Q4)-RF*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(&
        Q6)*COS(Q4))^2)/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4&
        ))^2)^0.5)*(U4*SIN(Q4)*(L1+L2*SIN(THETA)+L3*COS(THETA)*COS(Q6&
        )+2*RF*COS(THETA)*SIN(Q6)*(SIN(Q4)*COS(Q6)+2*SIN(THETA)*SIN(Q&
        6)*COS(Q4))/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2&
        )^0.5+RF*COS(THETA)*SIN(Q6)*(SIN(Q4)*COS(Q6)+2*SIN(THETA)*SIN&
        (Q6)*COS(Q4))*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4)&
        )^2)/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^1.5)+&
        COS(THETA)*COS(Q4)*(L3*U6*SIN(Q6)+RF*(SIN(Q6)*(SIN(Q4)*COS(Q6&
        )+SIN(THETA)*SIN(Q6)*COS(Q4))^2*(U6*SIN(Q4)*SIN(Q6)+SIN(THETA&
        )*U4*SIN(Q4)*SIN(Q6)-U4*COS(Q4)*COS(Q6)-SIN(THETA)*U6*COS(Q4)&
        *COS(Q6))-2*SIN(Q6)*(U4*COS(Q4)*COS(Q6)-U6*SIN(Q4)*SIN(Q6))*(&
        1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^2-2*U6*COS(&
        Q6)*(SIN(Q4)*COS(Q6)+2*SIN(THETA)*SIN(Q6)*COS(Q4))*(1-(SIN(Q4&
        )*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^2-SIN(Q6)*(SIN(Q4)*C&
        OS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^4*(U6*SIN(Q4)*SIN(Q6)+SIN(&
        THETA)*U4*SIN(Q4)*SIN(Q6)-U4*COS(Q4)*COS(Q6)-SIN(THETA)*U6*CO&
        S(Q4)*COS(Q6))-SIN(Q6)*(U4*COS(Q4)*COS(Q6)-U6*SIN(Q4)*SIN(Q6)&
        )*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)*(1-(SIN&
        (Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)-U6*COS(Q6)*(SIN(Q&
        4)*COS(Q6)+2*SIN(THETA)*SIN(Q6)*COS(Q4))*(-1+(SIN(Q4)*COS(Q6)&
        +SIN(THETA)*SIN(Q6)*COS(Q4))^2)*(1-(SIN(Q4)*COS(Q6)+SIN(THETA&
        )*SIN(Q6)*COS(Q4))^2))/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)&
        *COS(Q4))^2)^2.5))/(COS(Q4)^2*(L1+L2*SIN(THETA)+L3*COS(THETA)&
        *COS(Q6)+2*RF*COS(THETA)*SIN(Q6)*(SIN(Q4)*COS(Q6)+SIN(THETA)*&
        SIN(Q6)*COS(Q4))/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q&
        4))^2)^0.5+RF*COS(THETA)*SIN(Q6)*(SIN(Q4)*COS(Q6)+SIN(THETA)*&
        SIN(Q6)*COS(Q4))*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(&
        Q4))^2)/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^1.&
        5)^2) - (L2*COS(THETA)*U4*SIN(Q4)-RR*U4*SIN(Q4)-L3*(U4*SIN(Q6&
        )*COS(Q4)+U6*SIN(Q4)*COS(Q6)+SIN(THETA)*U4*SIN(Q4)*COS(Q6)+SI&
        N(THETA)*U6*SIN(Q6)*COS(Q4))-2*RF*(SIN(Q4)*COS(Q6)+SIN(THETA)&
        *SIN(Q6)*COS(Q4))*(U6*SIN(Q4)*SIN(Q6)+SIN(THETA)*U4*SIN(Q4)*S&
        IN(Q6)-U4*COS(Q4)*COS(Q6)-SIN(THETA)*U6*COS(Q4)*COS(Q6))/(1-(&
        SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^0.5-RF*(SIN(Q4&
        )*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))*(-1+(SIN(Q4)*COS(Q6)+SI&
        N(THETA)*SIN(Q6)*COS(Q4))^2)*(U6*SIN(Q4)*SIN(Q6)+SIN(THETA)*U&
        4*SIN(Q4)*SIN(Q6)-U4*COS(Q4)*COS(Q6)-SIN(THETA)*U6*COS(Q4)*CO&
        S(Q6))/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^1.5&
        )/(COS(Q4)*(L1+L2*SIN(THETA)+L3*COS(THETA)*COS(Q6)+2*RF*COS(T&
        HETA)*SIN(Q6)*(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))/(1&
        -(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^0.5+RF*COS(T&
        HETA)*SIN(Q6)*(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))*(-&
        1+(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)/(1-(SIN(Q4)&
        *COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^1.5)))*D2>
W_E_D>=0>
W_F_E>=U6*F3>
W_G_F>=U7*G2>
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
% BUILD THE POSITION VECTOR FROM CO TO DO
P_CO_DO>=L4*D1>+L5*D3>
% BUILD THE POSITION VECTOR FROM GO TO FO
P_GO_FO>=L6*F1>+L7*F3>
%---------------------------------------------------------------------%
%         DEFINE THE GENERALIZED SPEEDS
%---------------------------------------------------------------------%
Q1'=U1
Q2'=U2
Q3'=U3
Q4'=U4
Q5'=U5
Q6'=U6
Q7'=U7
%---------------------------------------------------------------------%
%         VELOCITIES
%---------------------------------------------------------------------%
% CALCULATE THE VELOCITY OF POINT CO
V_CO_N>=DT(P_NO_CO>,N)
% CALCULATE THE VELOCITY OF POINT CN
V2PTS(N,C,CO,CN)
% CALCULATE THE VELOCITY OF POINT GO
V_GO_N>=DT(P_NO_GO>,N)
% CALCULATE THE VELOCITY OF POINT GN
V2PTS(N,G,GO,GN)
% CALCULATE THE VELOCITY OF POINT DO
V2PTS(N,D,CO,DO)
% CALCULATE THE VELOCITY OF POINT FO
V2PTS(N,F,GO,FO)
%---------------------------------------------------------------------%
%         MOTION CONSTRAINTS
%---------------------------------------------------------------------%
% DUE TO THE ASSUMPTIONS OF NO SIDE SLIP AND NO SLIP ROLLING THE
% VELOCITIES OF THE FRONT AND REAR WHEEL CONTACT POINTS, CN AND GN,
% CANNOT HAVE COMPONENTS OF VELOCITY IN THE GROUND PLANE
DEPENDENT[1]=DOT(V_CN_N>,N1>)
DEPENDENT[2]=DOT(V_CN_N>,N2>)
DEPENDENT[3]=DOT(V_GN_N>,N1>)
DEPENDENT[4]=DOT(V_GN_N>,N2>)
% THE REAR WHEEL ANGULAR SPEED, U5, THE ROLL RATE, U4, AND THE
% STEERING RATE, U6, ARE TAKEN TO BE THE INDEPENDENT GENERALIZED
% SPEEDS
CONSTRAIN(DEPENDENT[U1,U2,U3,U7])
%---------------------------------------------------------------------%
%         ANGULAR ACCELERATIONS
%---------------------------------------------------------------------%
ALF_A_N>=DT(W_A_N>,N)
ALF_B_A>=DT(W_B_A>,N)
ALF_C_B>=DT(W_C_B>,N)
ALF_D_B>=DT(W_D_B>,N)
ALF_E_D>=DT(W_E_D>,N)
ALF_F_E>=DT(W_F_E>,N)
%---------------------------------------------------------------------%
%         ACCELERATIONS
%---------------------------------------------------------------------%
A_CO_N>=DT(V_CO_N>,N)
A_GO_N>=DT(V_GO_N>,N)
A2PTS(N,D,CO,DO)
A2PTS(N,F,GO,FO)
%---------------------------------------------------------------------%
%         FORCES
%---------------------------------------------------------------------%
GRAVITY(-G*N3>,C,D,F,G)
%---------------------------------------------------------------------%
%         EQUATIONS OF MOTION
%---------------------------------------------------------------------%
ZERO=FR()+FRSTAR()
KANE()
%---------------------------------------------------------------------%
%         SIMULATION
%---------------------------------------------------------------------%
UnitSystem  kg,meter,sec
INPUT U5=16.7 RAD/S
OUTPUT Q1 M,Q2 M,Q3 RAD,Q4 RAD,Q5 RAD,Q6 RAD,Q7 RAD
OUTPUT U1 M/S,U2 M/S,U3 RAD/S,U4 RAD/S,U5 RAD/S,U6 RAD/S,U7 RAD/S
CODE DYNAMICS() BICYCLE13.M
%---------------------------------------------------------------------%
CONSTANTS NU5 % NOMINAL SOLUTION FOR U5
%---------------------------------------------------------------------%
SOLVE(ZERO,U4',U5',U6')
Q4P=TAYLOR(Q4',1,Q1=0,Q2=0,Q3=0,Q4=0,Q5=0,Q6=0,Q7=0,U4=0,U5=NU5,U6=0)
Q6P=TAYLOR(Q6',1,Q1=0,Q2=0,Q3=0,Q4=0,Q5=0,Q6=0,Q7=0,U4=0,U5=NU5,U6=0)
PAUSE
A[1,1]=COEF(Q4P,Q4)
A[1,2]=COEF(Q4P,Q6)
A[1,3]=COEF(Q4P,U4)
A[1,4]=COEF(Q4P,U6)
A[2,1]=COEF(Q6P,Q4)
A[2,2]=COEF(Q6P,Q6)
A[2,3]=COEF(Q6P,U4)
A[2,4]=COEF(Q6P,U6)
%---------------------------------------------------------------------%
A[3,1]=EVALUATE(D(U4',Q4),Q1=0,Q2=0,Q3=0,Q4=0,Q5=0,Q6=0,Q7=0,U4=0,U5=NU5,U6=0)
PAUSE
A[3,2]=EVALUATE(D(U4',Q6),Q1=0,Q2=0,Q3=0,Q4=0,Q5=0,Q6=0,Q7=0,U4=0,U5=NU5,U6=0)
PAUSE
A[3,3]=EVALUATE(D(U4',U4),Q1=0,Q2=0,Q3=0,Q4=0,Q5=0,Q6=0,Q7=0,U4=0,U5=NU5,U6=0)
PAUSE
A[3,4]=EVALUATE(D(U4',U6),Q1=0,Q2=0,Q3=0,Q4=0,Q5=0,Q6=0,Q7=0,U4=0,U5=NU5,U6=0)
PAUSE
A[4,1]=EVALUATE(D(U6',Q4),Q1=0,Q2=0,Q3=0,Q4=0,Q5=0,Q6=0,Q7=0,U4=0,U5=NU5,U6=0)
PAUSE
A[4,2]=EVALUATE(D(U6',Q6),Q1=0,Q2=0,Q3=0,Q4=0,Q5=0,Q6=0,Q7=0,U4=0,U5=NU5,U6=0)
PAUSE
A[4,3]=EVALUATE(D(U6',U4),Q1=0,Q2=0,Q3=0,Q4=0,Q5=0,Q6=0,Q7=0,U4=0,U5=NU5,U6=0)
PAUSE
A[4,4]=EVALUATE(D(U6',U6),Q1=0,Q2=0,Q3=0,Q4=0,Q5=0,Q6=0,Q7=0,U4=0,U5=NU5,U6=0)
%---------------------------------------------------------------------%
%         SAVE OUTPUT
%---------------------------------------------------------------------%
SAVE BICYCLE13.ALL