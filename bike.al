% bike.al
AutoZ On
Newtonian N
UNITSYSTEM KG,METER,SECOND
BODIES A,B,C,D
VARIABLES Q{7}',qu8
CONSTANTS lambda, R, l1, h, a, S1A, S3A, S1B, S3B, MA, MB, MC, MD, g
MASS A=MA, B=MB, C=MC, D=MD
INERTIA A,IA11,IA22,IA33,IA12,IA23,IA31
INERTIA B,IB11,IB22,IB33,IB12,IB23,IB31
INERTIA C,IC11,IC22,IC11
INERTIA D,ID11,ID22,ID11
POINTS O, P1, P2, PSA
%
%
%
% ... GEOMETRY ...
%
% orientation of A in N
FRAMES A1,A2
SIMPROT(N,A1,3,Q3)
SIMPROT(A1,A2,1,Q4)
% Q8 is dependent on Q5 by the nonlinear geometric constraint dot(P_P1_P2>,n3>) = 0 which constraints the front
% wheel to stay on the ground. It is useful that Q8 does not appear in the Autolev code, therefore the equation
% Q8=f(Q4,Q5) is already inserted at this point
% dot(P_P1_P2>,n3>) = 0 linearized in Q8 around Q8=0 and solved for Q8:
qu8 =  -((h+a*sin(lambda))*cos(Q4)-h*cos(Q4)-R*cos(Q4)-a*(sin(Q4)*sin(Q5)+sin(lambda)*cos(Q4)*cos(Q5))-R*(-1+(sin(Q4)*cos(Q5)-sin(lambda)*sin(Q5)&
*cos(Q4))^2)/(1-(sin(Q4)*cos(Q5)-sin(lambda)*sin(Q5)*cos(Q4))^2)^0.5)/(cos(Q4)*(l1+sin(lambda)*(a*tan(lambda)+h/cos(lambda))+a*cos(lambda)*cos(Q5)-2*R*cos(lambda)*sin(Q5)&
*(sin(Q4)*cos(Q5)-sin(lambda)*sin(Q5)*cos(Q4))/(1-(sin(Q4)*cos(Q5)-sin(lambda)*sin(Q5)*cos(Q4))^2)^0.5-R*cos(lambda)*sin(Q5)*(sin(Q4)*cos(Q5)-sin(lambda)*sin(Q5)*cos(Q4))*&
(-1+(sin(Q4)*cos(Q5)-sin(lambda)*sin(Q5)*cos(Q4))^2)/(1-(sin(Q4)*cos(Q5)-sin(lambda)*sin(Q5)*cos(Q4))^2)^1.5))
%substitute for Q8 in equations above
SIMPROT(A2,A,2,qu8)
SIMPROT(A,C,2,Q6) % orientation of C in A
% orientation of B in A
FRAMES B1
SIMPROT(A,B1,2,-lambda) 
SIMPROT(B1,B,3,Q5)
SIMPROT(B,D,2,Q7) % orientation of D in B
P_O_P1>=Q1*N1>+Q2*N2> % position rear wheel contact point P1 in N
P_P1_CO> = R*A23> %position rear wheel axis point CO from rear wheel contact point P1
P_CO_PSA> = l1*A1>+h*A3> %position steering axis point PSA from rear rear wheel axis point CO
P_PSA_DO> = a*B1>-(h/cos(lambda)+a*tan(lambda))*B3> %position front wheel axis point DO from steering axis point PSA 
%find contact point P2
rf_direction> = UNITVEC(n3> - DOT(n3>,b2>)*b2>)
P_DO_P2> = -R*rf_direction> %position front wheel contact point P2 to front wheel axis point DO
P_CO_AO>=S1A*A1> + S3A*A3>
P_DO_BO>=S1B*A1> + S3B*A3>
%
%
%
% ... VELOCITIES ...
%
MotionVariables' U{7}'
Q1' = cos(Q3)*U1-cos(Q3)*R*sin(Q4)*U3-U2*sin(Q3)-R*cos(Q4)*U4*sin(Q3)
Q2' = cos(Q3)*U2+cos(Q3)*R*cos(Q4)*U4+sin(Q3)*U1-sin(Q3)*R*sin(Q4)*U3
Q3' = U3
Q4' = U4
Q5' = U5
Q6' = U6/R
Q7' = U7/R
W_A1_N> = U3*A13>
W_A2_A1> = U4*A21>
k1 = -(h*SIN(Q4)+R*SIN(Q4)-(h+a*SIN(lambda))*SIN(Q4)-a*(SIN(Q5)*COS(Q4)-SIN(lambda)*SIN(Q4)*COS(Q5))-2*R*(COS(Q4)*COS(Q5)+SIN(lambda)*SIN(Q4)*SIN(Q5))*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))/(1-(S&
IN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^0.5-R*(COS(Q4)*COS(Q5)+SIN(lambda)*SIN(Q4)*SIN(Q5))*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))*(-1+(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)/(&
1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^1.5)/(COS(Q4)*(l1+SIN(lambda)*(a*TAN(lambda)+h/COS(lambda))+a*COS(lambda)*COS(Q5)-2*R*COS(lambda)*SIN(Q5)*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(&
Q4))/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^0.5-R*COS(lambda)*SIN(Q5)*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))*(-1+(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)/(1-(SIN(Q4)*COS(Q&
5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^1.5))-((h+a*SIN(lambda))*COS(Q4)-h*COS(Q4)-R*COS(Q4)-a*(SIN(Q4)*SIN(Q5)+SIN(lambda)*COS(Q4)*COS(Q5))-R*(-1+(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)/(1-(SI&
N(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^0.5)*(SIN(Q4)*(l1+SIN(lambda)*(a*TAN(lambda)+h/COS(lambda))+a*COS(lambda)*COS(Q5)-2*R*COS(lambda)*SIN(Q5)*(SIN(Q4)*COS(Q5)-2*SIN(lambda)*SIN(Q5)*COS(Q4)&
)/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^0.5-R*COS(lambda)*SIN(Q5)*(SIN(Q4)*COS(Q5)-2*SIN(lambda)*SIN(Q5)*COS(Q4))*(-1+(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)/(1-(SIN(Q4)*COS(Q5&
)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^1.5)-R*COS(lambda)*SIN(Q5)*COS(Q4)*((COS(Q4)*COS(Q5)+SIN(lambda)*SIN(Q4)*SIN(Q5))*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^4-2*COS(Q4)*COS(Q5)*(1-(SIN(Q4)*COS(Q&
5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^2-(COS(Q4)*COS(Q5)+SIN(lambda)*SIN(Q4)*SIN(Q5))*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2-COS(Q4)*COS(Q5)*(-1+(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)&
*(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2))/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^2.5)/(COS(Q4)^2*(l1+SIN(lambda)*(a*TAN(lambda)+h/COS(lambda))+a*COS(lambda)*COS(Q5)-2*R*COS(&
lambda)*SIN(Q5)*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^0.5-R*COS(lambda)*SIN(Q5)*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))*(-1+(SIN(Q4)*&
COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^1.5)^2)
k2 = (SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))*(a-2*R*(SIN(Q4)*SIN(Q5)+SIN(lambda)*COS(Q4)*COS(Q5))/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^0.5-R*(SIN(Q4)*SIN(Q5)+SIN(lambda)*COS(Q4)*COS(Q&
5))*(-1+(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^1.5)/(COS(Q4)*(l1+SIN(lambda)*(a*TAN(lambda)+h/COS(lambda))+a*COS(lambda)*COS(Q5)-2*R*COS&
(lambda)*SIN(Q5)*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^0.5-R*COS(lambda)*SIN(Q5)*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))*(-1+(SIN(Q4)&
*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^1.5))-COS(lambda)*((h+a*SIN(lambda))*COS(Q4)-h*COS(Q4)-R*COS(Q4)-a*(SIN(Q4)*SIN(Q5)+SIN(lambda)*COS(Q4)*&
COS(Q5))-R*(-1+(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^0.5)*(a*SIN(Q5)+R*(2*COS(Q5)*(SIN(Q4)*COS(Q5)-2*SIN(lambda)*SIN(Q5)*COS(Q4))*(1-(SI&
N(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^2+COS(Q5)*(SIN(Q4)*COS(Q5)-2*SIN(lambda)*SIN(Q5)*COS(Q4))*(-1+(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)*(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*CO&
S(Q4))^2)+SIN(Q5)*((SIN(Q4)*SIN(Q5)+SIN(lambda)*COS(Q4)*COS(Q5))*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^4-2*SIN(Q4)*SIN(Q5)*(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^2-(SIN(Q4)*SIN(Q&
5)+SIN(lambda)*COS(Q4)*COS(Q5))*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2-SIN(Q4)*SIN(Q5)*(-1+(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)*(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2&
)))/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^2.5)/(COS(Q4)*(l1+SIN(lambda)*(a*TAN(lambda)+h/COS(lambda))+a*COS(lambda)*COS(Q5)-2*R*COS(lambda)*SIN(Q5)*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)&
*COS(Q4))/(1-(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^0.5-R*COS(lambda)*SIN(Q5)*(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))*(-1+(SIN(Q4)*COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)/(1-(SIN(Q4)*&
COS(Q5)-SIN(lambda)*SIN(Q5)*COS(Q4))^2)^1.5)^2)
W_A_A2> = (Z12*U4 + Z13*U5)*A2>
W_A_N> = W_A1_N> + W_A2_A1> + W_A_A2>
W_B_A> = U5*B3>
V_CO_N> = U1*A11> + U2*A12> - R*SIN(Q4)*U4*A13> %V_CO_N> = EXPRESS(Dt(P_O_CO>,N),A1)
W_C_A> = U6/R*A2>
W_D_B> = U7/R*B2>
V2PTS(N,A,CO,PSA)
V2PTS(N,B,PSA,DO)
V2PTS(N,A,CO,AO)
V2PTS(N,B,DO,BO)
%
%
%
% ... MOTION CONSTRAINTS ...
%
V2PTS(N,C,CO,P1)
V2PTS(N,D,DO,P2)
EXPRESS(V_P2_N>,N)
dependent[1] = DOT(V_P1_N>,N1>)
dependent[2] = DOT(V_P1_N>,N2>)
dependent[3] = DOT(V_P2_N>,N1>)
dependent[4] = DOT(V_P2_N>,N2>)
CONSTRAIN(DEPENDENT[U2,U3,U6,U7])
%
%
%
% ... FORCES ...
%
GRAVITY(-G*N3>)
%
%
%
% ... EQUATIONS OF MOTION ...
A_AO_N> = DT(V_AO_N>,N)
A_BO_N> = DT(V_BO_N>,N)
A_CO_N> = DT(V_CO_N>,N)
A_DO_N> = DT(V_DO_N>,N)
ZERO=FR()+FRSTAR()
KANE()
INPUT Q1=0 m, Q2=0 m, Q3=0 rad, Q4=0 rad, Q5=0 rad, Q6=0 rad, Q7=0 rad
INPUT U1=3 m/s, U4=0 rad/s, U5=0 rad/s
INPUT lambda=15/180*pi rad, R=0.5 m, l1 = 1.5 m, h=0.6 m, a=0.1 m, S1A=0.5 m, S3A=0.6 m, S1B=0.02 m, S3B=0.3 m, MA=80 kg, MB=2 kg, MC=2 kg, MD=2 kg
INPUT g=9.81 m/s^2
INPUT IA11=7 KG*M^2,IA22=8 KG*M^2,IA33=8 KG*M^2,IA12=0 KG*M^2,IA23=0 KG*M^2,IA31=0 KG*M^2
INPUT IB11=0.1 KG*M^2,IB22=0.1 KG*M^2,IB33=0 KG*M^2,IB12=0 KG*M^2,IB23=0 KG*M^2,IB31=0 KG*M^2
INPUT IC11=0.2 KG*M^2,IC22=0.5 KG*M^2
INPUT ID11=0.2 KG*M^2,ID22=0.5 KG*M^2
OUTPUT Q1, Q2, Q3, Q4, Q5, Q6, Q7, U1, U2 m/s, U3 rad/s, U4, U5, U6 m/s, U7 m/s, qu8 rad
%CODE DYNAMICS() model.m
%save out.all