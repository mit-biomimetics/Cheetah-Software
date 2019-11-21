close all
clc
clear all
%%
load('./../matlab_log/sim_data.mat')
fin = length(simulator_state.time);


%% Figure 1

% Formats
ylab = {'x (m)','y (m)','z (m)','Roll (deg)','Pitch (deg)','Yaw (deg)'};
ymax = [0.15 0.06 0.25 10 20 30];
ymin = -ymax;ymin(3) = 0.2;

figure(1)
for i = 1:3
   ax(i) = subplot(2,3,i);
   plot(simulator_state.time(:),CONTROLLER_two_contact_stand_data.p_act(1:fin,i))
   hold on
   plot(simulator_state.time(:),CONTROLLER_two_contact_stand_data.p_des(1:fin,i))
   xlabel('Time (s)')
   ylabel(ylab{i})
   axis([-inf inf ymin(i) ymax(i)])
   
   ax(i+3) = subplot(2,3,i+3);
   plot(simulator_state.time(:),rad2deg(CONTROLLER_two_contact_stand_data.rpy_act(1:fin,i)))
   hold on
   plot(simulator_state.time(:),rad2deg(CONTROLLER_two_contact_stand_data.rpy(1:fin,i)))
   xlabel('Time (s)')
   ylabel(ylab{i+3})
   axis([-inf inf ymin(i+3) ymax(i+3)])
end
linkaxes(ax,'x')

%% Figure 2

% Formats
ylab = {'F_x (m)','F_y (m)','F_z (m)'};
ymax = [10 10 60];
ymin = -ymax;ymin(3) = 0;
leg = {'FR','FL','BR','BL'};

figure(2)
for i = 1:4
   axf(i) = subplot(3,4,i);
   plot(simulator_state.time(:),CONTROLLER_two_contact_stand_data.f_control(1:fin,3*(i-1)+1))
   xlabel('Time (s)')
   ylabel(ylab{1})
   title(leg{i})
   axis([-inf inf ymin(1) ymax(1)])
   
   axf(i+4) = subplot(3,4,i+4);
   plot(simulator_state.time(:),CONTROLLER_two_contact_stand_data.f_control(1:fin,3*(i-1)+2))
   xlabel('Time (s)')
   ylabel(ylab{2})
   axis([-inf inf ymin(2) ymax(2)])
   
   axf(i+8) = subplot(3,4,i+8);
   plot(simulator_state.time(:),CONTROLLER_two_contact_stand_data.f_control(1:fin,3*(i-1)+3))
   xlabel('Time (s)')
   ylabel(ylab{3})
   axis([-inf inf ymin(3) ymax(3)])
end
linkaxes(axf,'x')

%% Figure 3 - comprhensive view

% Formats
ylab = {'Roll (deg)','Pitch (deg)','Yaw (deg)','x (m)','y (m)','z (m)','F_x (m)','F_y (m)','F_z (m)'};
ymax = [10 20 30 0.15 0.06 0.25 10 10 60];
ymin = [-10 -20 -30 -0.15 -0.06 0.2 -10 -10 0];

figure(3)
for i = 1:3
   axc(i) = subplot(3,3,i);
   plot(simulator_state.time(:),rad2deg(CONTROLLER_two_contact_stand_data.rpy_act(1:fin,i)))
   hold on
   plot(simulator_state.time(:),rad2deg(CONTROLLER_two_contact_stand_data.rpy(1:fin,i)))
   xlabel('Time (s)')
   ylabel(ylab{i})
   axis([-inf inf ymin(i) ymax(i)])
   
   axc(i+3) = subplot(3,3,i+3);
   plot(simulator_state.time(:),CONTROLLER_two_contact_stand_data.p_act(1:fin,i))
   hold on
   plot(simulator_state.time(:),CONTROLLER_two_contact_stand_data.p_des(1:fin,i))
   xlabel('Time (s)')
   ylabel(ylab{i+3})
   axis([-inf inf ymin(i+3) ymax(i+3)])
   
   axc(i+6) = subplot(3,3,i+6);
   plot(simulator_state.time(:),CONTROLLER_two_contact_stand_data.f_control(1:fin,i),'r')
   hold on
   plot(simulator_state.time(:),CONTROLLER_two_contact_stand_data.f_control(1:fin,i+3),'g--')
   plot(simulator_state.time(:),CONTROLLER_two_contact_stand_data.f_control(1:fin,i+6),'b--')
   plot(simulator_state.time(:),CONTROLLER_two_contact_stand_data.f_control(1:fin,i+9),'k')
   xlabel('Time (s)')
   ylabel(ylab{i+6})
   legend({'FR','FL','BR','BL'},'Location','NorthWest');
   axis([-inf inf ymin(i+6) ymax(i+6)])

end
linkaxes(axc,'x')




