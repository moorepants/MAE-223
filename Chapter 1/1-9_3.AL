% 1-9_3.AL
AUTOZ OFF
FRAMES A,B,C,D
POINTS O,P
VARIABLES Q{3}
CONSTANTS R
SIMPROT(A,B,-A2>,Q1)
SIMPROT(B,C,-B1>,Q2)
SIMPROT(C,D,-C2>,Q3)
P_O_P>=R*D1>
VARIABLES Q{3}'
IMPLICIT>=DT(P_O_P>,C)
Q3=3*T
Q3'=DT(Q3)
EXPLICIT>=EXPLICIT(EXPRESS(IMPLICIT>,C))
