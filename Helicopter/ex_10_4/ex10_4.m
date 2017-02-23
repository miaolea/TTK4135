
%A1,B1,u,x1,x2,x3,x4 = optCalc();

init;

h = 0.25;
alfa = 0.2;
beta = 20;
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

ul 	    = -30*pi/180;                   % Lower bound on control -- u1
uu 	    = 30*pi/180;                    % Upper bound on control -- u1

xl      = -Inf*ones(nx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(nx,1);               % Upper bound on states (no bound)
xl(3)   = ul;                           % Lower bound on state x3
xu(3)   = uu;                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = genbegr2(N, N, xl, xu, ul, uu); % hint: genbegr2
vlb(N*nx+N*nu)  = 0;                    % We want the last input to be zero
vub(N*nx+N*nu)  = 0;                    % We want the last input to be zero

Q1 = zeros(nx,nx);
Q1(1,1) = 1;
P1 = zeros(nu,nu); 
P1(1,1) = q_1;
P1(2,2) = q_2;
Q = 2*genq2(Q1,P1,N,N,nu);

f =@(z) 0.5*z'*Q*z;

%Aeq = gena2(A_d,B_d,N,nx,nu);
Aeq = [eye(N*nx) + kron(diag(ones(N-1,1),-1),-A),kron(eye(N),-B_d)];
%beq = zeros(1,size(Aeq,1));
beq = [-A_d*x0;zeros(nx*(N-1),1)];
%beq(1:nx) = -A_d*x0;

opts = optimset('MaxIter',15000, 'MaxFunEvals',20000);

tic
[z,fval] = fmincon(f,zeros(N*(nx+nu),1),[],[],Aeq,beq,vlb,vub,@constraint,opts); 
t1=toc;


phi1 = 0.0;
PhiOut = zeros(N*nx+N*nu,1);
for i=1:N*nx+N*nu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

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


