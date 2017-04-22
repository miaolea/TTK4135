%% Plotting

figure(2)
subplot(511)
stairs(t,simin(:,2)),grid
title('Estimates and measurements','FontSize',16);
ylabel('u')
subplot(512)
plot(t,x1,'m',t,x1,'mo',state_measurement2.time,state_measurement2.signals.values(:,1),'b'),grid
ylabel('lambda')
subplot(513)
plot(t,x2,'m',t,x2','mo',state_measurement2.time,state_measurement2.signals.values(:,2),'b'),grid
ylabel('r')
subplot(514)
plot(t,x3,'m',t,x3,'mo',state_measurement2.time,state_measurement2.signals.values(:,3),'b'),grid
ylabel('p')
subplot(515)
plot(t,x4,'m',t,x4','mo',state_measurement2.time,state_measurement2.signals.values(:,4),'b'),grid
xlabel('tid (s)'),ylabel('pdot')


%% Export to CSV
data = [state_measurement2.time 
        simin(:,2) 
        state_measurement2.signals.values(:,1) 
        state_measurement2.signals.values(:,2) 
        state_measurement2.signals.values(:,3) 
        state_measurement2.signals.values(:,4)];
csvwrite('ex10_2_measured_LQ.csv', data)