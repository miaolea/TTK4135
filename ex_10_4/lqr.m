Q = eye(6);

R = eye(2);

K = dlqr(A_d,B_d,Q,R,[]);

delta_t	= 0.25;
t = 0:delta_t:delta_t*(length(pc)-1);

uin = [t' pc ec];

xin = [t' x1 x2 x3 x4 x5 x6];