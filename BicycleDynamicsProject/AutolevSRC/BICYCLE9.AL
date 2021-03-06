% JASON MOORE -- MAY 22, 2006
AUTOZ ON
AUTORHS OFF
OVERWRITE ALL
% DECLARE THE INERTIAL REFERENCE FRAME
NEWTONIAN N
% DECLARE THREE INTERMEDIATE FRAMES
FRAMES A,B,E
% DECLARE FOUR BODIES: REAR WHEEL,FRAME,FORK,FRONT WHEEL
BODIES C,D,F,G
% DECLARE FOUR POINTS: FRONT CONTACT POINT,REAR CONTACT POINT (ONE ON GROUND AND ONE ON WHEEL)
POINTS NC,CN,NG,GN
% DECLARE THE FRAME DIMENSIONS
CONSTANTS RF,RR,L{7},THETA
% DECLARE THE GENERALIZED COORDINATES
VARIABLES Q{7}',PITCH
MOTIONVARIABLES' U{7}'

MASS C=MC,D=MD,F=MF,G=MG
INERTIA C,C1,C2,C3
INERTIA D,D1,D2,D3
INERTIA F,F1,F2,F3
INERTIA G,G1,G2,G3

% FRAME YAW
SIMPROT(N,A,3,Q3)
% FRAME ROLL
SIMPROT(A,B,-1,Q4)
% REAR WHEEL ROTATION
SIMPROT(B,C,2,Q5)
% FRAME PITCH
PITCH=(L2*COS(THETA)*COS(Q4)+L3*(SIN(Q4)*SIN(Q6)-SIN(THETA)*COS(Q4)*COS(&
        Q6))-RR*COS(Q4)-RF*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)/&
        (1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^0.5)/(COS(Q4)*(L1+&
        L2*SIN(THETA)+L3*COS(THETA)*COS(Q6)+2*RF*COS(THETA)*SIN(Q6)*(SIN(Q4)*&
        COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(&
        Q6)*COS(Q4))^2)^0.5+RF*COS(THETA)*SIN(Q6)*(SIN(Q4)*COS(Q6)+SIN(THETA)*&
        SIN(Q6)*COS(Q4))*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)/(1&
        -(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^1.5))
SIMPROT(B,D,-2,PITCH)
% HEAD TUBE ANGLE
SIMPROT(D,E,-2,THETA)
% STEERING ANGLE
SIMPROT(E,F,3,Q6)
% FRONT WHEEL ROTATION
SIMPROT(F,G,2,Q7)

% CALCULATE ANGULAR VELOCITIES
W_A_N>=U3*A3>
W_B_A>=-U4*B1>
W_C_B>=U5/RR*C2>
W_D_B>=-((L2*COS(THETA)*COS(Q4)+L3*(SIN(Q4)*SIN(Q6)-SIN(THETA)*COS(Q4)*COS&
        (Q6))-RR*COS(Q4)-RF*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)&
        /(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^0.5)*(U4*SIN(Q4)*(&
        L1+L2*SIN(THETA)+L3*COS(THETA)*COS(Q6)+2*RF*COS(THETA)*SIN(Q6)*(SIN(Q4)&
        *COS(Q6)+2*SIN(THETA)*SIN(Q6)*COS(Q4))/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*&
        SIN(Q6)*COS(Q4))^2)^0.5+RF*COS(THETA)*SIN(Q6)*(SIN(Q4)*COS(Q6)+2*SIN(&
        THETA)*SIN(Q6)*COS(Q4))*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4)&
        )^2)/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^1.5)+COS(THETA)&
        *COS(Q4)*(L3*U6*SIN(Q6)+RF*(SIN(Q6)*(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)&
        *COS(Q4))^2*(U6*SIN(Q4)*SIN(Q6)+SIN(THETA)*U4*SIN(Q4)*SIN(Q6)-U4*COS(&
        Q4)*COS(Q6)-SIN(THETA)*U6*COS(Q4)*COS(Q6))-2*SIN(Q6)*(U4*COS(Q4)*COS(&
        Q6)-U6*SIN(Q4)*SIN(Q6))*(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))&
        ^2)^2-2*U6*COS(Q6)*(SIN(Q4)*COS(Q6)+2*SIN(THETA)*SIN(Q6)*COS(Q4))*(1-(&
        SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^2-SIN(Q6)*(SIN(Q4)*COS(&
        Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^4*(U6*SIN(Q4)*SIN(Q6)+SIN(THETA)*U4*SIN&
        (Q4)*SIN(Q6)-U4*COS(Q4)*COS(Q6)-SIN(THETA)*U6*COS(Q4)*COS(Q6))-SIN(Q6)*&
        (U4*COS(Q4)*COS(Q6)-U6*SIN(Q4)*SIN(Q6))*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)&
        *SIN(Q6)*COS(Q4))^2)*(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)&
        -U6*COS(Q6)*(SIN(Q4)*COS(Q6)+2*SIN(THETA)*SIN(Q6)*COS(Q4))*(-1+(SIN(Q4)&
        *COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)*(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*&
        SIN(Q6)*COS(Q4))^2))/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)&
        ^2.5))/(COS(Q4)^2*(L1+L2*SIN(THETA)+L3*COS(THETA)*COS(Q6)+2*RF*COS(THE&
        TA)*SIN(Q6)*(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))/(1-(SIN(Q4)*&
        COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^0.5+RF*COS(THETA)*SIN(Q6)*(SIN(&
        Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)&
        *SIN(Q6)*COS(Q4))^2)/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)&
        ^1.5)^2) - (L2*COS(THETA)*U4*SIN(Q4)-RR*U4*SIN(Q4)-L3*(U4*SIN(Q6)*COS(&
        Q4)+U6*SIN(Q4)*COS(Q6)+SIN(THETA)*U4*SIN(Q4)*COS(Q6)+SIN(THETA)*U6*SIN(&
        Q6)*COS(Q4))-2*RF*(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))*(U6*SIN(&
        Q4)*SIN(Q6)+SIN(THETA)*U4*SIN(Q4)*SIN(Q6)-U4*COS(Q4)*COS(Q6)-SIN(THETA)&
        *U6*COS(Q4)*COS(Q6))/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)&
        ^0.5-RF*(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))*(-1+(SIN(Q4)*COS(&
        Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)*(U6*SIN(Q4)*SIN(Q6)+SIN(THETA)*U4*&
        SIN(Q4)*SIN(Q6)-U4*COS(Q4)*COS(Q6)-SIN(THETA)*U6*COS(Q4)*COS(Q6))/(1-(&
        SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^1.5)/(COS(Q4)*(L1+L2*SIN&
        (THETA)+L3*COS(THETA)*COS(Q6)+2*RF*COS(THETA)*SIN(Q6)*(SIN(Q4)*COS(Q6)+&
        SIN(THETA)*SIN(Q6)*COS(Q4))/(1-(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(&
        Q4))^2)^0.5+RF*COS(THETA)*SIN(Q6)*(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*&
        COS(Q4))*(-1+(SIN(Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)/(1-(SIN(&
        Q4)*COS(Q6)+SIN(THETA)*SIN(Q6)*COS(Q4))^2)^1.5)))*D2>
W_E_D>=0>
W_F_E>=U6*F3>
W_G_F>=U7/RF*G2>

% CALCULATE THE UNIT VECTOR ALONG THE LINE FROM GO TO GN
RFHAT>=UNITVEC(DOT(F2>,N3>)*F2>-N3>)
% BUILD THE POSITION VECTOR FROM THE INERTIAL ORIGIN TO NC
P_NO_NC>=Q1*N1>+Q2*N2>
% BUILD THE POSITION VECTOR FROM NC TO CN
P_NC_CN>=0>
% BUILD THE POSITION VECTOR FROM CN TO CO
P_CN_CO>=RR*B3>
% BUILD THE POSITION VECTOR FROM CO TO GO
P_CO_GO>=L1*D1>-L2*F3>+L3*F1>
% BUILD THE POSITION VECTOR FROM GO TO GN
P_GO_GN>=RF*RFHAT>
% BUILD THE POSITION VECTOR FROM GN TO NG
P_GN_NG>=0>
% BUILD THE POSITION VECTOR FROM CO TO DO
P_CO_DO>=L4*D1>+L5*D3>
% BUILD THE POSITION VECTOR FROM GO TO FO
P_GO_FO>=L6*F1>+L7*F3>

% DEFINE THE GENERALIZED SPEEDS
Q1'=COS(Q3)*(U1+RR*SIN(Q4)*U3)-SIN(Q3)*(U2-RR*COS(Q4)*U4)
Q2'=SIN(Q3)*(U1+RR*SIN(Q4)*U3)+COS(Q3)*(U2-RR*COS(Q4)*U4)
Q3'=U3
Q4'=U4
Q5'=U5/RR
Q6'=U6
Q7'=U7/RF

% CALCULATE THE VELOCITY OF POINT CO
V_CO_N>=-RR*U3*SIN(Q4)*B1>+RR*U4*B2>+(COS(Q3)*(U1+RR*U3*SIN(Q4))-SIN(Q3)*(U2-RR*U4*COS(Q4)))*N1> + U2*N2>
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

%DEP1=DOT(V_CN_N>,N1>)
%DEP2=DOT(V_CN_N>,N2>)
%DEP3=DOT(V_GN_N>,N1>)
%DEP4=DOT(V_GN_N>,N2>)
%SOLVE([DEP1;DEP2;DEP3;DEP4],U2,U3,U5,U7)

%DEPENDENT[1]=DOT(V_CN_N>,N1>)
%DEPENDENT[2]=DOT(V_CN_N>,N2>)
%DEPENDENT[3]=DOT(V_GN_N>,N1>)
%DEPENDENT[4]=DOT(V_GN_N>,N2>)

DEPENDENT[1]=EXPLICIT(DOT(V_CN_N>,N1>))
DEPENDENT[2]=EXPLICIT(DOT(V_CN_N>,N2>))
DEPENDENT[3]=EXPLICIT(DOT(V_GN_N>,N1>))
DEPENDENT[4]=EXPLICIT(DOT(V_GN_N>,N2>))

%DEPENDENT=EXPLICIT(DEPENDENT,Q1',Q2',Q3',Q4',Q5',Q6',Q7')

%DEPENDENT[3]=REPLACE(DOT(V_GN_N>,N1>),Q1'=(COS(Q3)*(U1+RR*SIN(Q4)*U3)-SIN(Q3)*(U2-RR*COS(Q4)*U4)))
%DEPENDENT[4]=REPLACE(DOT(V_GN_N>,N2>),Q2'=(SIN(Q3)*(U1+RR*SIN(Q4)*U3)+COS(Q3)*(U2-RR*COS(Q4)*U4)))

CONSTRAIN(DEPENDENT[U2,U3,U5,U7])

%ALF_A_N>=DT(W_A_N>,N)
%ALF_B_A>=DT(W_B_A>,N)
%ALF_C_B>=DT(W_C_B>,N)
%ALF_D_B>=DT(W_D_B>,N)
%ALF_E_D>=DT(W_E_D>,N)
%ALF_F_E>=DT(W_F_E>,N)

%GRAVITY(-G*N3>,C,D,F,G)

% SAVE OUTPUT
SAVE BICYCLE9.ALL