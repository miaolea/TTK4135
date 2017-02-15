
%A1,B1,u,x1,x2,x3,x4 = optCalc();

Q = eye(4);
R = 0.1;

K = dlqr(A1,B1,Q,R,[]);

delta_t	= 0.25;
t = 0:delta_t:delta_t*(length(u)-1);

uin = [t' u];

xin = [t' x1 x2 x3 x4];