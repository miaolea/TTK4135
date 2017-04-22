
%% Initialize system
init;

h = 0.25;
q_1 = 1;
q_2 = 1;
N = 40;
nx = 6;
nu = 2;


A = [0 1 0 0 0 0;
    0 0 -K_2 0 0 0;
    0 0 0 1 0 0;
    0 0 -K_1*K_pp -K_1*K_pd 0 0;
    0 0 0 0 0 1;
    0 0 0 0 -K_3*K_ep -K_3*K_ed];

B = [0 0;
    0 0;
    0 0;
    K_1*K_pp 0;
    0 0;
    0 K_3*K_ep];

A_d = h.*A + eye(6);
B_d = h.*B;

x1_0 = pi;
x2_0 = 0;
x3_0 = 0;
x4_0 = 0;
x5_0 = 0;
x6_0 = 0;
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';

z = zeros(N*(nx+nu),1);

%% Generate constraints on measurements and inputs
ul 	    = [-40*pi/180; -Inf];                   
uu 	    = [40*pi/180; Inf];                   

xl      = -Inf*ones(nx,1);            
xu      = Inf*ones(nx,1);             

% xl(2)   = [-45*pi/180];
% xu(2)   = [45*pi/180];
% xl(6)   = [-7*pi/180];
% xu(6)   = [7*pi/180];
xl(3)   = ul(1);                         
xu(3)   = uu(1);                          


[vlb,vub]       = genbegr2(N, N, xl, xu, ul, uu);
vlb(N*nx+N*nu)  = 0;                  
vub(N*nx+N*nu)  = 0;                   

%% Generate objective funtion and constraints
Q = kron(eye(N),diag([1 zeros(1,nx-1)]));
R = kron(eye(N),diag([q_1 q_2]));

G = blkdiag(Q,R);

f =@(z) 0.5*z'*G*z;

Aeq1 = eye(N*nx) -kron(diag(diag(eye(N-1)),-1), A_d);
Aeq2 = kron(diag(diag(eye(N))),-B_d);
Aeq = [Aeq1 Aeq2];

beq = [A_d*x0;zeros(nx*(N-1),1)];

%% Calculate optimal trajectory and input
opts = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',15000);

tic
[z,fval] = fmincon(f,[x0(1); zeros(N*(nx+nu)-1,1)],[],[],Aeq,beq,vlb,vub,@constraint,opts); 
t1=toc;


%% Extract control inputs and states
pc  = [0;z(N*nx+1:nu:N*(nx+nu))];
ec  = [0;z(N*nx+2:nu:N*(nx+nu))];

x1 = [x0(1);z(1:nx:N*nx)];
x2 = [x0(2);z(2:nx:N*nx)];
x3 = [x0(3);z(3:nx:N*nx)];
x4 = [x0(4);z(4:nx:N*nx)];
x5 = [x0(5);z(5:nx:N*nx)];
x6 = [x0(6);z(6:nx:N*nx)];

num_variables = (nx+nu)/h;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

pc   = [zero_padding; pc; zero_padding];
ec   = [zero_padding; ec; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];

%% LQR
lQr = diag([2 1 1 2 10 1]);
lQr(5,5) = 10;

lqR = [1 0; 0 0.1];

K = dlqr(A_d,B_d,lQr,lqR,[]);

delta_t	= 0.25;
t = 0:delta_t:delta_t*(length(pc)-1);

uin = [t' pc ec];

xin = [t' x1 x2 x3 x4 x5 x6];

%csvwrite('ex10_4_optimal.csv', xin)

figure(1)
subplot(1,1,1)
hold on;
stairs(t,pc,'r'),grid
stairs(t,ec,'m'),grid
ylabel('u')


figure(2)
subplot(3,2,1)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
title('Estimates and measurements','FontSize',16);

subplot(3,2,2)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')

subplot(3,2,3)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')

subplot(3,2,4)
plot(t,x4,'m',t,x4','mo'),grid
ylabel('pdot')

subplot(3,2,5)
plot(t,x5,'m',t,x5','mo'),grid
xlabel('tid (s)'),ylabel('e')

subplot(3,2,6)
plot(t,x6,'m',t,x6','mo'),grid
xlabel('tid (s)'),ylabel('edot')