close all
clc
clear all
%%
addpath('./functions')
%load('./../matlab_log/data_exp_run5_freq.mat')
% load('./../matlab_log/hallway2.mat')
% load('./../matlab_log/exp_hallway.mat')
% load('./../matlab_log/exp_treadmill1.mat')
%load('./../matlab_log/data_exp_bounding.mat')
% load('./../matlab_log/data_exp_bounding.mat')
load('./../matlab_log/run_3.mat')

fig = fn_open_figures(6);

%%
st_idx = 4000;
end_idx = length(leg_control_data.lcm_timestamp);%-85000;
%end_idx = st_idx + 1000; %length(wbc_lcm_data.lcm_timestamp);%-85000;
%time = linspace(0, 1, length(wbc_lcm_data.lcm_timestamp));
time = leg_control_data.lcm_timestamp;

figure(fig(1))
for i =1:12
    subplot(4,3,i)
hold on
plot(time(st_idx:end_idx), leg_control_data.tau_est(st_idx:end_idx,i))
plot(time(st_idx:end_idx), leg_control_command.tau_ff(st_idx:end_idx,i))
grid on
axis tight
end
xlabel('Tau')

% JPos
figure(fig(2))
for i =1:12
    subplot(4,3,i)
hold on
plot(time(st_idx:end_idx), leg_control_data.q(st_idx:end_idx,i))
plot(time(st_idx:end_idx), leg_control_command.q_des(st_idx:end_idx,i))
grid on
axis tight
end
xlabel('JPos')

% JPos
figure(fig(3))
for i =1:12
    subplot(4,3,i)
hold on
plot(time(st_idx:end_idx), leg_control_data.qd(st_idx:end_idx,i))
plot(time(st_idx:end_idx), leg_control_command.qd_des(st_idx:end_idx,i))
grid on
axis tight
end
xlabel('JVel')


% foot pos
figure(fig(4))
for i =1:12
    subplot(4,3,i)
hold on
plot(time(st_idx:end_idx), leg_control_data.p(st_idx:end_idx,i))
grid on
axis tight
end
xlabel('Foot pos')

% foot vel
figure(fig(5))
for i =1:4
    subplot(4,1,i)
hold on
plot(leg_control_data.lcm_timestamp, leg_control_data.v(:,3*i-2))
plot(leg_control_data.lcm_timestamp, leg_control_data.v(:,3*i))
plot(wbc_lcm_data.lcm_timestamp, wbc_lcm_data.Fr(:,3*i)*0.05) 

grid on
axis tight
xlim([70, 80])
end
xlabel('Foot vel')


% 
figure(fig(6))
for i =1:12
    subplot(4,3,i)
plot(time(st_idx:end_idx), leg_control_command.kp_joint(st_idx:end_idx,i))
grid on
axis tight
end
xlabel('gains')

