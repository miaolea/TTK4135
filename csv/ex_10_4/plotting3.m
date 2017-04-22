% Plotting
t = 0:h:h*(length(pc)-1);

% figure(2)
% subplot(7,1,1)
% hold on;
% stairs(t,pc,'r'),grid
% stairs(t,ec,'m'),grid
% title('q1 = 1, q2 = 1','FontSize',16);
% ylabel('u')
% subplot(7,1,2)
% plot(t,x1,'m',t,x1,'mo'),grid
% ylabel('lambda')
% subplot(7,1,3)
% plot(t,x2,'m',t,x2','mo'),grid
% ylabel('r')
% subplot(7,1,4)
% plot(t,x3,'m',t,x3,'mo'),grid
% ylabel('p')
% subplot(7,1,5)
% plot(t,x4,'m',t,x4','mo'),grid
% xlabel('tid (s)'),ylabel('pdot')
% subplot(7,1,6)
% plot(t,x5,'m',t,x5','mo'),grid
% xlabel('tid (s)'),ylabel('e')
% subplot(7,1,7)
% plot(t,x6,'m',t,x6','mo'),grid
% xlabel('tid (s)'),ylabel('edot')

% figure(2)
% subplot(2,1,1)
% hold on;
% stairs(t,pc,'r'),grid
% stairs(t,ec,'m'),grid
% title('Optimal and actual input','FontSize',16);
% ylabel('u optimal')
% 
% subplot(2,1,2)
% hold on;
% stairs(input_actual.time, input_actual.signals.values(:,1),'b'),grid
% stairs(input_actual.time, input_actual.signals.values(:,2),'c'),grid
% ylabel('u actual')
% 
% print -depsc ex10_4_input.eps

figure(1)
subplot(1,1,1)
hold on;
stairs(t,pc,'r'),grid
stairs(t,ec,'m'),grid
ylabel('u')

%% Export to CSV
data1 = [t' pc ec];
csvwrite('ex10_4_u_fb_elevconst_t_pc_ec.csv', data1)

figure(2)
subplot(3,2,1)
plot(t,x1,'m',t,x1,'mo',state_measurement3_fb_travel_const.time, state_measurement3_fb_travel_const.signals.values(:,1),'b'),grid
ylabel('lambda')
title('Estimates and measurements','FontSize',16);

subplot(3,2,2)
plot(t,x2,'m',t,x2','mo',state_measurement3_fb_travel_const.time,state_measurement3_fb_travel_const.signals.values(:,2),'b'),grid
ylabel('r')

subplot(3,2,3)
plot(t,x3,'m',t,x3,'mo',state_measurement3_fb_travel_const.time,state_measurement3_fb_travel_const.signals.values(:,3),'b'),grid
ylabel('p')

subplot(3,2,4)
plot(t,x4,'m',t,x4','mo',state_measurement3_fb_travel_const.time,state_measurement3_fb_travel_const.signals.values(:,4),'b'),grid
ylabel('pdot')

subplot(3,2,5)
plot(t,x5,'m',t,x5','mo',state_measurement3_fb_travel_const.time,state_measurement3_fb_travel_const.signals.values(:,5),'b'),grid
xlabel('tid (s)'),ylabel('e')

subplot(3,2,6)
plot(t,x6,'m',t,x6','mo',state_measurement3_fb_travel_const.time,state_measurement3_fb_travel_const.signals.values(:,6),'b'),grid
xlabel('tid (s)'),ylabel('edot')

%% Export to CSV
data2 = [state_measurement3_fb_travel_const.time 
        state_measurement3_fb_travel_const.signals.values(:,1) 
        state_measurement3_fb_travel_const.signals.values(:,2) 
        state_measurement3_fb_travel_const.signals.values(:,3) 
        state_measurement3_fb_travel_const.signals.values(:,4)
        state_measurement3_fb_travel_const.signals.values(:,5)
        state_measurement3_fb_travel_const.signals.values(:,6)];
%csvwrite('ex10_4_measured_fb_elevconst_t_states.csv', data2)


