
Q = diag([10 1 1 1]);

R = 1;

K = dlqr(A_d,B_d,Q,R,[]);

delta_t	= 0.25;
t = 0:delta_t:delta_t*(length(u)-1);

uin = [t' u];

xin = [t' x1 x2 x3 x4];