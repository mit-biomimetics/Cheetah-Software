clear all
clc

rpy = [0.3, 0.1, 0.5];
Rx = [1 0 0; 
    0 cos(rpy(1)) -sin(rpy(1)); 
    0 sin(rpy(1)) cos(rpy(1))];
Ry = [cos(rpy(2)) 0 sin(rpy(2)); 
    0 1 0; 
    -sin(rpy(2)) 0 cos(rpy(2))];
Rz = [cos(rpy(3)) -sin(rpy(3)) 0; 
    sin(rpy(3)) cos(rpy(3)) 0; 
    0 0 1];

Rot = Rz*Ry*Rx;

shoulder = [0.25, 0.17, 0].';
body_pos = [1.2, 0.3, 0.3].';
body_vel = [0.2, 0.1, 0.05].';
body_ang_vel = [0.1, 0.3, 0.2].';
step_time = 0.25;

foot_loc = body_pos + Rot*shoulder + step_time/2 * (body_vel + cross(body_ang_vel, Rot*shoulder));
foot_loc = foot_loc + body_pos(2)/9.81*cross(body_vel, body_ang_vel);
foot_loc