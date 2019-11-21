close all
clc
clear all
%%
load('./../matlab_log/sim_data.mat')
fin = length(CONTROLLER_two_contact_stand_data.p_act);


%% Figure 1

% Formats
ylab = {'x (m)','y (m)','z (m)','Roll (deg)','Pitch (deg)','Yaw (deg)'};
ymax = [0.15 0.06 0.25 10 20 30];
ymin = -ymax;ymin(3) = 0.2;

figure
for i = 1:3
    ax(i) = subplot(2,3,i);
    plot(CONTROLLER_two_contact_stand_data.p_act(1:fin,i))
    hold on
    plot(CONTROLLER_two_contact_stand_data.p_des(1:fin,i))
    xlabel('Time (s)')
    ylabel(ylab{i})
    axis([-inf inf ymin(i) ymax(i)])
    
    ax(i+3) = subplot(2,3,i+3);
    plot(rad2deg(CONTROLLER_two_contact_stand_data.rpy_act(1:fin,i)))
    hold on
    plot(rad2deg(CONTROLLER_two_contact_stand_data.rpy(1:fin,i)))
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

% Find the point in f where the contact state for OOC feet dips below
% threshold
i = 1;
threshold = 0.1;
trans = length(CONTROLLER_two_contact_stand_data.contact_state(:,1));
while i < length(CONTROLLER_two_contact_stand_data.contact_state(:,1))
    if CONTROLLER_two_contact_stand_data.contact_state(i,1) < threshold
        trans = i;
        break
    end
    i = i+1;
end

figure
for i = 1:4
    axf(i) = subplot(3,4,i);
    hold on
    plot(CONTROLLER_two_contact_stand_data.f_control(1:fin,3*(i-1)+1),'r')
    plot(CONTROLLER_two_contact_stand_data.f_control(1:trans-1,3*(i-1)+1),'b')
    xlabel('Time (s)')
    ylabel(ylab{1})
    title(leg{i})
    axis([-inf inf ymin(1) ymax(1)])
    
    axf(i+4) = subplot(3,4,i+4);
    hold on
    plot(CONTROLLER_two_contact_stand_data.f_control(1:fin,3*(i-1)+2),'r')
    plot(CONTROLLER_two_contact_stand_data.f_control(1:trans-1,3*(i-1)+2),'b')
    xlabel('Time (s)')
    ylabel(ylab{2})
    axis([-inf inf ymin(2) ymax(2)])
    
    axf(i+8) = subplot(3,4,i+8);
    hold on
    plot(CONTROLLER_two_contact_stand_data.f_control(1:fin,3*(i-1)+3),'r')
    plot(CONTROLLER_two_contact_stand_data.f_control(1:trans-1,3*(i-1)+3),'b')
    plot(CONTROLLER_two_contact_stand_data.maxForces(1:fin,i),'k--')
    plot(CONTROLLER_two_contact_stand_data.f_ref(1:fin,i),'g:')
    xlabel('Time (s)')
    ylabel(ylab{3})
    axis([-inf inf ymin(3) ymax(3)])
end
linkaxes(axf,'x')

figure
for i = 1:4
    subplot(1,5,i)
    plot(CONTROLLER_two_contact_stand_data.contact_state(:,i))
end
subplot(1,5,5)
plot(CONTROLLER_two_contact_stand_data.stance_legs);


%% Figure 3 - comprhensive view

% Formats
ylab = {'Roll (deg)','Pitch (deg)','Yaw (deg)','x (m)','y (m)','z (m)','F_x (m)','F_y (m)','F_z (m)'};
ymax = [10 20 30 0.15 0.06 0.25 10 10 60];
ymin = [-10 -20 -30 -0.15 -0.06 0.2 -10 -10 0];

figure
for i = 1:3
    axc(i) = subplot(3,3,i);
    plot(rad2deg(CONTROLLER_two_contact_stand_data.rpy_act(1:fin,i)))
    hold on
    plot(rad2deg(CONTROLLER_two_contact_stand_data.rpy(1:fin,i)))
    xlabel('Time (s)')
    ylabel(ylab{i})
    axis([-inf inf ymin(i) ymax(i)])
    
    axc(i+3) = subplot(3,3,i+3);
    plot(CONTROLLER_two_contact_stand_data.p_act(1:fin,i))
    hold on
    plot(CONTROLLER_two_contact_stand_data.p_des(1:fin,i))
    xlabel('Time (s)')
    ylabel(ylab{i+3})
    axis([-inf inf ymin(i+3) ymax(i+3)])
    
    axc(i+6) = subplot(3,3,i+6);
    plot(CONTROLLER_two_contact_stand_data.f_control(1:fin,i),'r')
    hold on
    plot(CONTROLLER_two_contact_stand_data.f_control(1:fin,i+3),'g--')
    plot(CONTROLLER_two_contact_stand_data.f_control(1:fin,i+6),'b--')
    plot(CONTROLLER_two_contact_stand_data.f_control(1:fin,i+9),'k')
    xlabel('Time (s)')
    ylabel(ylab{i+6})
    legend({'FR','FL','BR','BL'},'Location','NorthWest');
    axis([-inf inf ymin(i+6) ymax(i+6)])
    
end
linkaxes(axc,'x')




