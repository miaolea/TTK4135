
%% Initialize system
init;
h	= 0.25;

A = [0 1 0 0 ; 0 0 -K_2 0 ; 0 0 0 1 ; 0 0 -K_1*K_pp -K_1*K_pd];
B = [0 ; 0 ; 0 ; K_1*K_pp];

nx = size(A,2);
nu = size(B,2);

A_d = eye(nx)+h*A;
B_d = h*B;

x1_0 = pi;
x2_0 = 0;
x3_0 = 0;
x4_0 = 0;
x0 = [x1_0 x2_0 x3_0 x4_0]';

q = 0.1;
N  = 60;                             
z  = zeros(N*(nx+nu),1);             
z0 = z;                            

%% Generate constraints on measurements and inputs
ul 	    = -30*pi/180;
uu 	    = 30*pi/180;

xl      = -Inf*ones(nx,1);
xu      = Inf*ones(nx,1);
xl(3)   = ul;
xu(3)   = uu;

[vlb,vub]       = genbegr2(N, N, xl, xu, ul, uu);
vlb(N*nx+N*nu)  = 0;
vub(N*nx+N*nu)  = 0;

%% Generate the matrix Q and the vector c
%% (objecitve function weights in the QP problem) 
Q = kron(eye(N),diag([1 zeros(1,nx-1)]));
R = kron(eye(N),q);
c = zeros(N*(nx+nu),1);
G = 2*blkdiag(Q,R);

%% Generate system matrixes for linear model
Aeq1 = eye(N*nx) - kron(diag(diag(eye(N-1)),-1),A_d);
Aeq2 = kron(diag(diag(eye(N))),-B_d);
Aeq = [Aeq1 Aeq2];

beq = zeros(N*nx,1);
beq(1:nx) = A_d*x0;

%% Solve QP
tic
[z,lambda] = quadprog(G,c,[],[],Aeq,beq,vlb,vub,x0);
t1=toc;

%% Extract control inputs and states
u  = [z(N*nx+1:N*nx+N*nu);z(N*nx+N*nu)]; 

x1 = [x0(1);z(1:nx:N*nx)];
x2 = [x0(2);z(2:nx:N*nx)];
x3 = [x0(3);z(3:nx:N*nx)];
x4 = [x0(4);z(4:nx:N*nx)];

num_variables = 5/h;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);


u   = [zero_padding; u; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];


%% Plotting
t = 0:h:h*(length(u)-1);

figure(1)
subplot(511)
stairs(t,u),grid
title('q = 0.1','FontSize',18);
ylabel('u')
subplot(512)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(513)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(514)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(515)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')

%% Export as CSV

data = [t' u x1 x2 x3 x4];
csvwrite('ex10_2_q01.csv', data)


%% Generate the matrix Q and the vector c
%% (objecitve function weights in the QP problem)
q = 10;
R = kron(eye(N),q);
G = 2*blkdiag(Q,R);


%% Solve QP
tic;
[z,lambda] = quadprog(G,c,[],[],Aeq,beq,vlb,vub,x0);
t1=toc;

%% Extract control inputs and states
u  = [z(N*nx+1:N*nx+N*nu);z(N*nx+N*nu)]; 

x1 = [x0(1);z(1:nx:N*nx)];
x2 = [x0(2);z(2:nx:N*nx)];
x3 = [x0(3);z(3:nx:N*nx)];
x4 = [x0(4);z(4:nx:N*nx)];

num_variables = 5/h;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);


u   = [zero_padding; u; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];


%% Plotting
t = 0:h:h*(length(u)-1);

figure(2)
subplot(511)
stairs(t,u),grid
title('q = 10','FontSize',18);
ylabel('u')
subplot(512)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(513)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(514)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(515)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')


%% Export as CSV

data = [t' u x1 x2 x3 x4];
csvwrite('ex10_2_q10.csv', data)



%% Generate the matrix Q and the vector c
%% (objecitve function weights in the QP problem) 
q = 1;
Q = kron(eye(N),diag([1 zeros(1,nx-1)]));
R = kron(eye(N),q);
c = zeros(N*(nx+nu),1);
G = 2*blkdiag(Q,R);

%% Generate system matrixes for linear model
Aeq1 = eye(N*nx) - kron(diag(diag(eye(N-1)),-1),A_d);
Aeq2 = kron(diag(diag(eye(N))),-B_d);
Aeq = [Aeq1 Aeq2];

beq = zeros(N*nx,1);
beq(1:nx) = A_d*x0;

%% Solve QP
tic
[z,lambda] = quadprog(G,c,[],[],Aeq,beq,vlb,vub,x0);
t1=toc;

%% Extract control inputs and states
u  = [z(N*nx+1:N*nx+N*nu);z(N*nx+N*nu)]; 

x1 = [x0(1);z(1:nx:N*nx)];
x2 = [x0(2);z(2:nx:N*nx)];
x3 = [x0(3);z(3:nx:N*nx)];
x4 = [x0(4);z(4:nx:N*nx)];

num_variables = 5/h;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);


u   = [zero_padding; u; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];

%% Plotting
t = 0:h:h*(length(u)-1);

figure(3)
subplot(511)
stairs(t,u),grid
title('q = 1','FontSize',18);
ylabel('u')
subplot(512)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(513)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(514)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(515)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')

%% Export as CSV

data = [t' u x1 x2 x3 x4];
csvwrite('ex10_2_q1.csv', data)

%%Export input
simin = [t' u];