%% Converts the GRF to Momentum rate of change
function hDot = NonlinearInput(u, r, s)

% Number of possible robot contact feet
NUM_FEET = 2;

% Set up symbolic CoM forces and torque vectors
f_com = sym(zeros(2, 1));
tau_com = sym(zeros(1, 1));

% Iterate through all of the robot legs
for foot = 1:NUM_FEET
    
    % Leg state number
    n = size(u,1)/NUM_FEET*(foot - 1);
    
    % Foot position vector from CoM (CURRENTLY ONLY FOR FLAT GROUND)
    r_foot = [r(1 + n);1;r(2 + n)];
    
    % Foot ground reaction forces
    f_foot = [u(1 + n);1;u(2 + n)];
    
    % Summation of CoM forces
    f_com = f_com + s(foot)*[f_foot(1);f_foot(3)];
    
    % Summation of CoM torques (take only tau_yy)
    tau_cross = CrossProd(r_foot)*f_foot;
    tau_com = tau_com + s(foot)*tau_cross(2);
end

% Net wrench on the CoM
hDot = [f_com; tau_com];

end