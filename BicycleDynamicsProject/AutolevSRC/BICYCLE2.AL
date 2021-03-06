NEWTONIAN N
FRAMES A,B,E
BODIES C,D,F,G
POINTS O,PF,PR
CONSTANTS RF,RR,L{3},THETA,D01,D03,FO1,FO3
VARIABLES Q{5}',Q7',Q8'
Q6 = (L1*COS(Q4)+RR*COS(Q4)+L2*SIN(THETA)*COS(Q4)-L3*(SIN(Q4)*SIN(Q7)-COS(THETA)*COS(Q4)*COS(Q7))-RR*COS(ASIN(SIN(Q4)*COS(Q7)+COS(THETA)*SIN(Q7)*COS(Q4))))/(COS(Q4)*(L2*COS(THETA)-L3*SIN(THETA)*COS(Q7)-RR*SIN(THETA)*SIN(Q7)*(SIN(Q4)*COS(Q7)+COS(THETA)*SIN(Q7)*COS(Q4))/(1-(SIN(Q4)*COS(Q7)+COS(THETA)*SIN(Q7)*COS(Q4))^2)^0.5))
SIMPROT(N,A,3,Q3)
SIMPROT(A,B,-1,Q4)
SIMPROT(B,C,2,Q5)
SIMPROT(B,D,2,Q6)
SIMPROT(D,E,-2,THETA)
SIMPROT(E,F,-1,Q7)
SIMPROT(F,G,-2,Q8)
I3>=DOT(-N3>,F2>)*DOT(F2>,N1>)*COS(ASIN(DOT(-N3>,F2>)))/(1+DOT(-N3>,F2>)*DOT(F2>,N3>))*N1>+DOT(-N3>,F2>)*DOT(F2>,N2>)*COS(ASIN(DOT(-N3>,F2>)))/(1+DOT(-N3>,F2>)*DOT(F2>,N3>))*N2>+COS(ASIN(DOT(-N3>,F2>)))*N3>
I2>=F2>
I1>=CROSS(I2>,I3>)
P_O_PR>=Q1*N1>+Q2*N2>
P_PR_PF>=RR*B3>+L1*D3>+L2*F1>+L3*F3>-RR*I3>
P_PR_CO>=RR*B3>
P_PR_DO>=RR*B3>+D1O*D1>+D03*D3>
P_PR_FO>=RR*B3>+L1*D3>+FO1*F1>+FO3*F3>
