%% Plotting
t = 0:h:h*(length(pc)-1);

figure(2)
subplot(7,1,1)
hold on;
stairs(t,pc,'r'),grid
stairs(t,ec,'m'),grid
title('q1 = 1, q2 = 1','FontSize',16);
ylabel('u')
subplot(7,1,2)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(7,1,3)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(7,1,4)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(7,1,5)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')
subplot(7,1,6)
plot(t,x5,'m',t,x5','mo'),grid
xlabel('tid (s)'),ylabel('e')
subplot(7,1,7)
plot(t,x6,'m',t,x6','mo'),grid
xlabel('tid (s)'),ylabel('edot')
