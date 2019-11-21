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
% load('./../matlab_log/sim_data.mat')
% load('./../matlab_log/0903_walking.mat')
% load('./../matlab_log/0903_body.mat')
load('./../matlab_log/0903_imu_1.mat')

fig = fn_open_figures(4);

%%
st_idx = 1;
end_idx = length(microstrain.lcm_timestamp);%-85000;
%end_idx = st_idx + 1000; %length(wbc_lcm_data.lcm_timestamp);%-85000;
%time = linspace(0, 1, length(wbc_lcm_data.lcm_timestamp));
time = microstrain.lcm_timestamp;
time_se = state_estimator.lcm_timestamp;
end_idx_se = length(time_se);

figure(fig(1))
for i =1:4
    subplot(4,1,i)
hold on
plot(time, microstrain.quat(st_idx:end_idx,i));
plot(time_se, state_estimator.quat(st_idx:end_idx_se,i));
grid on
axis tight
end
xlabel('Tau')

% JPos
figure(fig(2))
for i =1:3
    subplot(3,1,i)
hold on
plot(time, microstrain.omega(st_idx:end_idx, i))
plot(time_se, state_estimator.omegaBody(st_idx:end_idx_se,i))
grid on
axis tight
end
xlabel('JPos')


figure(fig(3))
for i =1:3
    subplot(3,1,i)
hold on
plot(time, microstrain.rpy(st_idx:end_idx,i));
plot(time_se, state_estimator.rpy(st_idx:end_idx_se,i));
grid on
axis tight
end
xlabel('Tau')



figure(fig(4))
for i =1:3
    subplot(3,1,i)
hold on
plot(time, microstrain.acc(st_idx:end_idx,i));
% plot(time, microstrain.bad_packets(st_idx:end_idx));
plot(hw_vectornav.lcm_timestamp, hw_vectornav.a(:,i));
grid on
axis tight
end
xlabel('Tau')

figure
plot(microstrain.bad_packets);
