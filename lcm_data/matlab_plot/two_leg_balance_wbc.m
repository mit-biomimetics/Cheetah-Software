clear;close all;clc

load('./../matlab_log/sim_data.mat')

fin = length((wbc_lcm_data.Fr_des(:,1)));


%% Foot Positions
% Formats
ylab = {'x (m)','y (m)','z (m)'};
col = {'r','b','g','k'};
figure
for i = 1:3
    ax(i) = subplot(1,3,i);
    hold on
    for j = [0 3 6 9]
    plot(wbc_lcm_data.foot_pos(:,i+j),col{j/3+1})
    plot(wbc_lcm_data.foot_pos_cmd(:,i+j),[col{j/3+1},'o'])
    end
    xlabel('Time (s)')
    ylabel(ylab{i})
    %axis([-inf inf ymin(i) ymax(i)])
end
linkaxes(ax,'x')

%% Forces
% Formats
ylab = {'F_x (m)','F_y (m)','F_z (m)'};
ymax = [10 10 60];
ymin = -ymax;ymin(3) = 0;
leg = {'FR','FL','BR','BL'};

figure
for i = 1:4
    axf(i) = subplot(3,4,i);
    hold on
    plot(wbc_lcm_data.Fr_des(1:fin,3*(i-1)+1),'--')
    plot(wbc_lcm_data.Fr(1:fin,3*(i-1)+1))
    xlabel('Time (s)')
    ylabel(ylab{1})
    title(leg{i})
    axis([-inf inf ymin(1) ymax(1)])
    
    axf(i+4) = subplot(3,4,i+4);
    hold on
    plot(wbc_lcm_data.Fr_des(1:fin,3*(i-1)+2),'--')
    plot(wbc_lcm_data.Fr(1:fin,3*(i-1)+2))
    xlabel('Time (s)')
    ylabel(ylab{2})
    axis([-inf inf ymin(2) ymax(2)])
    
    axf(i+8) = subplot(3,4,i+8);
    hold on
    plot(wbc_lcm_data.Fr_des(1:fin,3*(i-1)+3),'--')
    plot(wbc_lcm_data.Fr(1:fin,3*(i-1)+3))
    xlabel('Time (s)')
    ylabel(ylab{3})
    axis([-inf inf ymin(3) ymax(3)])
end
linkaxes(axf,'x')