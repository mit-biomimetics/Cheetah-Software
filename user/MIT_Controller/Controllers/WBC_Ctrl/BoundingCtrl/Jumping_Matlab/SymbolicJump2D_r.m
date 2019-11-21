%%
addpath Utils
clear; clc;
SymbolicJump2D_timer = tic;
params_timer = tic;
fprintf('General symbolic variable setup...\n')

%% Parameter setup
% Number of feet (front 2 and back 2 are paired together)
NUM_FEET = 2;

% Number of timesteps
N = 1;

% Decision Variables (for current timestep)
syms x z phi dx dz dphi Fxf Fzf Fxb Fzb rxf rxb real
states = [x;z;phi;dx;dz;dphi];
inputs = [Fxf;Fzf;rxf;Fxb;Fzb;rxb];

% Decision Variables (for previous timestep)
syms x0 z0 phi0 dx0 dz0 dphi0 Fxf0 Fzf0 Fxb0 Fzb0 rxf0 rxb0 real
states0 = [x0;z0;phi0;dx0;dz0;dphi0];
inputs0 = [Fxf0;Fzf0;rxf0;Fxb0;Fzb0;rxb0];

% Decision Variables (for next timestep)
syms x1 z1 phi1 dx1 dz1 dphi1 Fxf1 Fzf1 Fxb1 Fzb1 rxf1 rxb1 real
states1 = [x1;z1;phi1;dx1;dz1;dphi1];
inputs1 = [Fxf1;Fzf1;rxf1;Fxb1;Fzb1;rxb1];

% Desired states and inputs
syms x_d z_d phi_d dx_d dz_d dphi_d Fxf_d Fzf_d Fxb_d Fzb_d rxf_d rxb_d real
states_ref = [x_d;z_d;phi_d;dx_d;dz_d;dphi_d];
inputs_ref = [Fxf_d;Fzf_d;rxf_d;Fxb_d;Fzb_d;rxb_d];

% Decision Variable minimums
syms x_min z_min phi_min dx_min dz_min dphi_min Fxf_min Fzf_min Fxb_min Fzb_min rxf_min rxb_min real
states_min = [x_min;z_min;phi_min;dx_min;dz_min;dphi_min];
inputs_min = [Fxf_min;Fzf_min;rxf_min;Fxb_min;Fzb_min;rxb_min];

% Decision Variable maximums
syms x_max z_max phi_max dx_max dz_max dphi_max Fxf_max Fzf_max Fxb_max Fzb_max rxf_max rxb_max real
states_max = [x_max;z_max;phi_max;dx_max;dz_max;dphi_max];
inputs_max = [Fxf_max;Fzf_max;rxf_max;Fxb_max;Fzb_max;rxb_max];

% Footstep variables
syms pxf pzf pxb pzb rxf rzf rxb rzb real
r_foot = [rxf;rzf;rxb;rzb];
%p_foot = [pxf;pzf;pxb;pzb];

% Contact Variables
syms sf sb sf0 sb0 real
c_states = [sf;sb];
c_states0 = [sf0;sb0];

% Sparsity pattern parameters
syms iter NUM_X NUM_C real

% Optimization Parameters
NUM_STATES = size(states,1);
NUM_INPUTS = size(inputs,1);
NUM_DECISION_VARS = NUM_STATES + NUM_INPUTS;
NUM_CONSTRAINTS = NUM_STATES + 2*NUM_FEET + NUM_FEET;

% Time and environment parameters
syms dt g mu_g real

% Physical robot Parameters
syms m Iyy leg_length_max real
I_tensor = diag([m, m, Iyy]);
I_tensor_inv = I_tensor^-1;

% Robot linearized system
System.A = [eye(3), dt*eye(3); zeros(3), eye(3)];
System.B = [(dt^2)/2*I_tensor_inv; dt*I_tensor_inv];
System.C = eye(NUM_STATES);
System.D = zeros(NUM_STATES,NUM_INPUTS);
System.Dw = [((dt)^2)/2*eye(3); eye(3)*dt];
System.w = [0;g;0];

% Cost Function Weights
Q =  diag(sym('q%d', [NUM_STATES, 1], 'real'));  % states
R =  diag(sym('r%d', [NUM_INPUTS, 1], 'real'));  % inputs

% Lagrange multipliers on the hessian
lambda = sym('lambda%d', [NUM_CONSTRAINTS,1], 'real');

fprintf('DONE general symbolic variable setup: %f s\n\n',toc(params_timer));


%% Initialization
if true%INITIALIZATION
    % Start the timer for calculating constraints
    initialization_timer = tic;
    fprintf('Calculating Initialization...\n')
    
    % Set the infinity value
    INF = 2e19;
    
    for k = 1:N
        %%% STATE BOUNDS %%%
        % State upper bound
        states_ub = states_max;
        
        % State lower bound
        states_lb = states_min;
        
        for foot = 1:NUM_FEET
            % Leg state number
            n = NUM_INPUTS/NUM_FEET*(foot - 1);
            
            % Force input reference
            inputs_ref0(1+n:2+n,:) = 0*c_states(foot)*m*g;
            
            %%% INPUT BOUNDS %%%
            % Ground reaction force input upper bound
            inputs_ub(n+1:n+2,:) = c_states(foot)*inputs_max(n+1:n+2);
            
            % Ground reaction force input lower bound
            inputs_lb(n+1:n+2,:) = c_states(foot)*inputs_min(n+1:n+2);
            
        end
        %%% INITIAL STATES %%%
        % Convert GRF to momentum rate of change
        h_dot = NonlinearInput_r(states(:,k), inputs(:,k), c_states(:,k));
        
        % Compute simplified discrete dynamics
        y = simplify(DiscreteDynamics(states0(:,k), h_dot, System));
        
        %%% CONSTRAINT BOUNDS %%%
        % Constraint index
        i_c = NUM_CONSTRAINTS*(k - 1);
        
        % Dynamics, equality constraint {NUM_STATES}
        constraints_ub(i_c+1:i_c+NUM_STATES,1) = zeros(NUM_STATES,1);
        constraints_lb(i_c+1:i_c+NUM_STATES,1) = zeros(NUM_STATES,1);
        i_c = i_c + NUM_STATES;
        
        % Friction cone, inequality constraint {2*NUM_FEET}
        constraints_ub(i_c+1:i_c+2*NUM_FEET,1) = zeros(2*NUM_FEET,1);
        constraints_lb(i_c+1:i_c+2*NUM_FEET,1) = -INF*ones(2*NUM_FEET,1);
        i_c = i_c + 2*NUM_FEET;
    end
    
    % Initial decision variables
    decision_vars0 = [y; inputs_ref0];
    
    % Decision variable bounds
    decision_vars_lb = simplify([states_lb;inputs_lb]);
    decision_vars_ub = simplify([states_ub;inputs_ub]);
end


%% Cost, Gradient, and Cost Hessian
if true%COST
    % Start the timer for calculating cost
    cost_timer = tic;
    fprintf('Calculating Cost...\n')
    
    % Initialize cost
    cost = 0;
    for k = 1:N
        
        % Convert GRF to momentum rate of change
        hDot = NonlinearInput_r(states(:,k), inputs(:,k), c_states(:,k));
        
        % Compute simplified discrete dynamics
        y = DiscreteDynamics(states(:,k), hDot, System);
        
        % State trajectory tracking error
        x_error = (states_ref - states(:,k));
        %x_error = (states_ref - y);
        
        % Policy reference regularization
        u_error = (inputs_ref - inputs(:,k));
        
        % Objective cost function
        cost = cost + x_error'*Q*x_error + u_error'*R*u_error;
        
    end
    % Cost
    J = simplify(cost);
    
    % Decision variables
    decision_vars = [states;inputs];
    
    % Gradient
    fprintf('Calculating Cost Gradient...\n')
    cost_gradient = simplify(jacobian(J, decision_vars)');
    
    % Cost Hessian
    fprintf('Calculating Cost Hessian...\n')
    syms obj_factor real
    cost_hessian = simplify(obj_factor*jacobian(cost_gradient, decision_vars));
    
    % Diagonal entries of weight matrices only for generation
    Q =  diag(Q); R =  diag(R);
    
    % Print timing statistics
    fprintf('DONE calculating cost function: %f s\n\n',toc(cost_timer));
end


%% Constraints

%% Constraints
if true%CONSTRAINTS
    % Start the timer for calculating constraints
    constraint_timer = tic;
    fprintf('Calculating Constraints...\n')
    
    % Initialize iteration constraint vector
    constraints = sym(zeros(NUM_CONSTRAINTS, N));
    
    for k = 1:N
        % Constraint index
        i_c = NUM_CONSTRAINTS*(k - 1);
        
        % Convert GRF to momentum rate of change
        hDot = NonlinearInput_r(states(:,k), inputs(:,k), c_states(:,k));
        
        % Compute simplified discrete dynamics
        y = DiscreteDynamics(states0(:,k), hDot, System);
        
        % Dynamics, x_{k+1} - f(x_k, u_k) = 0 {NUM_STATES}
        constraints(i_c + 1:i_c + NUM_STATES,k) = states(:,k) - y;
        %constraints(i_c + 1:i_c + NUM_STATES,k) = states1(:,k) - y;
        i_c = i_c + NUM_STATES;
        
        % Add the constraints for each foot
        for foot = 1:NUM_FEET
            % Reset constraint index
            i_c = NUM_CONSTRAINTS*(k - 1) + NUM_STATES;
            
            % Leg state number
            n = NUM_INPUTS/NUM_FEET*(foot - 1);
            
            % Foot ground reaction forces
            f_foot = inputs(1+n:2+n, k);
            
            % Footstep vector previous
            rx_foot0 = inputs0(3+n, k);
            
            % Foot position in the world previous
            px_foot0 = (x0 + rx_foot0);
            
            % Footstep vector
            rx_foot = inputs(3+n, k);
            
            % Foot position in the world
            px_foot = (x + rx_foot);
            
            %% Input constraints
            % Friction pyramids {2*NUM_FEET}
            i_friction = i_c + 2*(foot - 1);
            constraints(i_friction + 1:i_friction + 2,k) = ...
                c_states(foot)*[...
                f_foot(2) - mu_g*f_foot(2);...
                -f_foot(2) - mu_g*f_foot(2)];
            i_c = i_c + 2*NUM_FEET;
            
            % No slip {NUM_FEET}
            i_slip = i_c + (foot - 1);
            constraints(i_slip + 1,k) = ...
                c_states0(foot)*c_states(foot)*(px_foot0 - px_foot);
            i_c = i_c + NUM_FEET;
            
        end
    end
    
    % Simplify the symbolic constraints
    constraints = simplify(constraints);
    
    % Calculate the constraint jacobian
    %    rows: constraints
    %    columns: decision variables
    fprintf('Calculating Constraint Jacobian...\n')
    constraint_jacobian = simplify(jacobian(constraints,[states0;inputs0;states;inputs]));
    
    % Non zero entries in the constraint jacobian
    constraint_jacobian_nz = constraint_jacobian(constraint_jacobian~=0);
    
    % Calculate the constraint jacobian sparsity pattern
    fprintf('Calculating Constraint Jacobian Sparsity Pattern...\n')
    [row_index_CJ, col_index_CJ] = find(constraint_jacobian);
    
    % Shift indices for C++
    row_index_CJ = row_index_CJ - 1;
    col_index_CJ = col_index_CJ - 1;
    
    % Add the iteration index modifier
    row_index_CJ = row_index_CJ + ones(size(row_index_CJ,1),1)*iter*NUM_C;
    col_index_CJ = col_index_CJ + ones(size(col_index_CJ,1),1)*iter*NUM_X;
    
    % Calculate the constraint hessian
    fprintf('Calculating Constraint Hessian...\n')
    constraint_hessian = zeros(2*NUM_DECISION_VARS:2*NUM_DECISION_VARS);
    
    % Sum of the hessian of each constraint with Lagrange multiplier
    for c = 1:NUM_CONSTRAINTS
        constraint_hessian = constraint_hessian + ...
            lambda(c)*jacobian(constraint_jacobian(c,:),[states0;inputs0;states;inputs]);
    end
    
    % Non zeros only in current timestep entries
    constraint_hessian = simplify(constraint_hessian(NUM_DECISION_VARS+1:2*NUM_DECISION_VARS,NUM_DECISION_VARS+1:2*NUM_DECISION_VARS));
    
    %% Initial Constraints
    % The first step does not include the previous dynamics as the initial
    % state is not a decision variable. Therefore it is a special case that
    % only relies on the current timestep
    fprintf('Calculating Initial Constraints...\n')
    
    % Copy the constraints
    constraints_initial = constraints;
    
    % Initial iteration constraint jacobian (no previous decision variables)
    fprintf('Calculating Final Constraint Jacobian...\n')
    constraint_jacobian_initial = simplify(jacobian(constraints_initial,[states;inputs]));
    
    % Non zero entries in the final constraint jacobian
    constraint_jacobian_initial_nz = ...
        constraint_jacobian_initial(constraint_jacobian_initial~=0);
    
    % Calculate the final constraint jacobian sparsity pattern
    fprintf('Calculating Initial Constraint Jacobian Sparsity Pattern...\n')
    [row_index_initial_CJ, col_index_initial_CJ] = ...
        find(constraint_jacobian_initial);
    
    % Shift indices for C++
    row_index_initial_CJ = row_index_initial_CJ - 1;
    col_index_initial_CJ = col_index_initial_CJ - 1;
    
    % Add the iteration index modifier
    row_index_initial_CJ = row_index_initial_CJ + ...
        ones(size(row_index_initial_CJ,1),1)*iter*NUM_C;
    col_index_initial_CJ = col_index_initial_CJ + ...
        ones(size(col_index_initial_CJ,1),1)*iter*NUM_X;
    
    % Calculate the initial constraint hessian
    fprintf('Calculating Initial Constraint Hessian...\n')
    constraint_hessian_initial = zeros(NUM_DECISION_VARS:NUM_DECISION_VARS);
    
    % Sum of the hessian of each constraint with Lagrange multiplier
    for c = 1:NUM_CONSTRAINTS
        constraint_hessian_initial = constraint_hessian_initial + ...
            lambda(c)*jacobian(constraint_jacobian_initial(c,:),[states;inputs]);
    end
    
    % Non zeros only in current timestep entries
    constraint_hessian_initial = ...
        simplify(constraint_hessian_initial(1:NUM_DECISION_VARS,1:NUM_DECISION_VARS));
    
    %% Final constraints
    % Since the final timestep does not include future dynamics or foot
    % locations, it needs to be treated as a special case where only the
    % instantaneous forces and foot vectors are constrained.
    fprintf('Calculating Final Constraints...\n')
    
    % Final constraints don't include future decision variables
    constraints_final = constraints;
    
    % Final dynamics are not constrained
    %constraints_final(1:NUM_STATES,1) = zeros(NUM_STATES,1);
    
    % Final iteration constraint jacobian (no future decision variables)
    fprintf('Calculating Final Constraint Jacobian...\n')
    constraint_jacobian_final = simplify(jacobian(constraints_final,[states;inputs]));
    
    % Non zero entries in the final constraint jacobian
    constraint_jacobian_final_nz = ...
        constraint_jacobian_final(constraint_jacobian_final~=0);
    
    % Calculate the final constraint jacobian sparsity pattern
    fprintf('Calculating Final Constraint Jacobian Sparsity Pattern...\n')
    [row_index_final_CJ, col_index_final_CJ] = ...
        find(constraint_jacobian_final);
    
    % Shift indices for C++
    row_index_final_CJ = row_index_final_CJ - 1;
    col_index_final_CJ = col_index_final_CJ - 1;
    
    % Add the iteration index modifier
    row_index_final_CJ = row_index_final_CJ + ...
        ones(size(row_index_final_CJ,1),1)*iter*NUM_C;
    col_index_final_CJ = col_index_final_CJ + ...
        ones(size(col_index_final_CJ,1),1)*iter*NUM_X;
    
    % Calculate the final constraint hessian
    fprintf('Calculating Final Constraint Hessian...\n')
    constraint_hessian_final = zeros(NUM_DECISION_VARS:NUM_DECISION_VARS);
    
    % Sum of the hessian of each constraint with Lagrange multiplier
    for c = 1:NUM_CONSTRAINTS
        constraint_hessian_final = constraint_hessian_final + ...
            lambda(c)*jacobian(constraint_jacobian_final(c,:),[states;inputs]);
    end
    
    % Non zeros only in current timestep entries
    constraint_hessian_final = ...
        simplify(constraint_hessian_final(1:NUM_DECISION_VARS,1:NUM_DECISION_VARS));
    
    % Print out the non zero entry results
    fprintf('   Nonzeros in Constraint Jacobian: %i\n',nnz(constraint_jacobian));
    fprintf('   Nonzeros in initial Constraint Jacobian: %i\n',nnz(constraint_jacobian_initial));
    fprintf('   Nonzeros in final Constraint Jacobian: %i\n',nnz(constraint_jacobian_final));
    
    % Print timing statistics
    fprintf('DONE calculating constraints: %f s\n\n',toc(constraint_timer));
    
end


%% Lagrangian Hessian
if true%LAGRANGIAN_HESSIAN
    hessian_timer = tic;
    fprintf('Calculating Lagrangian Hessian...\n');
    
    % Combine the cost and constraint hessians
    lagrangian_hessian = cost_hessian + constraint_hessian;
    
    % Only use the lower triangular part since it is symmetric
    lagrangian_hessian_full = lagrangian_hessian;
    lagrangian_hessian = simplify(tril(lagrangian_hessian));
    
    % Non zero entries in the lagrangian hessian
    lagrangian_hessian_nz = lagrangian_hessian(lagrangian_hessian~=0);
    
    fprintf('Calculating Lagrangian Hessian Sparsity Pattern...\n')
    
    % Find the non-zero row and column indices
    [row_index_H, col_index_H] = find(lagrangian_hessian);
    
    % Shift indices for C++
    row_index_H = row_index_H - 1;
    col_index_H = col_index_H - 1;
    
    % Add the iteration index modifier
    row_index_H = row_index_H + ones(size(row_index_H,1),1)*iter*NUM_X;
    col_index_H = col_index_H + ones(size(col_index_H,1),1)*iter*NUM_X;
    
    %% Initial Hessian
    fprintf('Calculating Initial Lagrangian Hessian...\n');
    
    % Combine the cost and constraint hessians
    lagrangian_hessian_initial = cost_hessian + constraint_hessian_initial;
    
    % Only use the lower triangular part since it is symmetric
    lagrangian_hessian_full_initial = lagrangian_hessian_initial;
    lagrangian_hessian_initial = simplify(tril(lagrangian_hessian_initial));
    
    % Non zero entries in the lagrangian hessian
    lagrangian_hessian_initial_nz = lagrangian_hessian_initial(lagrangian_hessian_initial~=0);
    
    fprintf('Calculating Initial Lagrangian Hessian Sparsity Pattern...\n')
    
    % Find the non-zero row and column indices
    [row_index_initial_H, col_index_initial_H] = find(lagrangian_hessian_initial);
    
    % Shift indices for C++
    row_index_initial_H = row_index_initial_H - 1;
    col_index_initial_H = col_index_initial_H - 1;
    
    % Add the iteration index modifier
    row_index_initial_H = row_index_initial_H +...
        ones(size(row_index_initial_H,1),1)*iter*NUM_X;
    col_index_initial_H = col_index_initial_H +...
        ones(size(col_index_initial_H,1),1)*iter*NUM_X;
    
    %% Final Hessian
    fprintf('Calculating Final Lagrangian Hessian...\n');
    
    % Combine the cost and constraint hessians
    lagrangian_hessian_final = cost_hessian + constraint_hessian_final;
    
    % Only use the lower triangular part since it is symmetric
    lagrangian_hessian_full_final = lagrangian_hessian_final;
    lagrangian_hessian_final = simplify(tril(lagrangian_hessian_final));
    
    % Non zero entries in the lagrangian hessian
    lagrangian_hessian_final_nz = lagrangian_hessian_final(lagrangian_hessian_final~=0);
    
    fprintf('Calculating Final Lagrangian Hessian Sparsity Pattern...\n')
    
    % Find the non-zero row and column indices
    [row_index_final_H, col_index_final_H] = find(lagrangian_hessian_final);
    
    % Shift indices for C++
    row_index_final_H = row_index_final_H - 1;
    col_index_final_H = col_index_final_H - 1;
    
    % Add the iteration index modifier
    row_index_final_H = row_index_final_H +...
        ones(size(row_index_final_H,1),1)*iter*NUM_X;
    col_index_final_H = col_index_final_H +...
        ones(size(col_index_final_H,1),1)*iter*NUM_X;
    
    % Print out the non zero entry results
    fprintf('   Nonzeros in Lagrangian Hessian: %i\n',nnz(lagrangian_hessian));
    fprintf('   Nonzeros in final Lagrangian Hessian: %i\n',nnz(lagrangian_hessian_final));
    
    fprintf('DONE calculating Lagrangian Hessian: %f s\n\n',toc(hessian_timer));
    
end



%% Generate Functions
if true%GENERATE_FUNCTIONS
    functions_timer = tic;
    fprintf('Generating MATLAB functions...\n')
    
    % Initialization
    
    % Generate decision variable and constraint bounds
    fprintf('Bounds function Gen...\n')
    matlabFunction(decision_vars_lb, decision_vars_ub, constraints_ub, constraints_lb,...
        'file','Jump2D/Jump2DBounds',...
        'vars',{states_max, states_min, inputs_max, inputs_min, c_states});
    
    % Generate initial guess, reference policy, and footstep references
    fprintf('Initial Conditions function Gen...\n')
    matlabFunction(decision_vars0, inputs_ref0,...
        'file','Jump2D/Jump2DInitialize',...
        'vars',{states, states0, inputs, c_states, dt, r_foot, m, Iyy, g});
    
    % Cost function
    
    % Generate cost function
    fprintf('Cost function Gen...\n');
    matlabFunction(J,...
        'file','Jump2D/Jump2DCost',...
        'vars',{states, inputs, states_ref, inputs_ref, Q, R});
    
    % Cost gradient
    
    % Generate cost gradient function
    fprintf('Cost Gradient function Gen...\n');
    matlabFunction(cost_gradient,...
        'file','Jump2D/Jump2DCostGradient',...
        'vars',{states, inputs, states_ref, inputs_ref, Q, R});
    
    % Constraints
    
    % Generate constraints function
    fprintf('Constraints function Gen...\n');
    matlabFunction(constraints,...
        'file','Jump2D/Jump2DConstraints',...
        'vars',{states, inputs, states0, inputs0, c_states, c_states0, dt,...
        r_foot, m, Iyy, g, mu_g});
    
    % Generate final constraints function
    fprintf('Final Constraints function Gen...\n');
    matlabFunction(constraints_final,...
        'file','Jump2D/Jump2DConstraintsFinal',...
        'vars',{states, inputs, states0, inputs0, c_states, c_states0, dt,...
        r_foot, m, Iyy, g, mu_g});
    
    % Generate initial constraints function
    fprintf('Initial Constraints function Gen...\n');
    matlabFunction(constraints_initial,...
        'file','Jump2D/Jump2DConstraintsInitial',...
        'vars',{states, inputs, states0, inputs0, c_states, c_states0, dt,...
        r_foot, m, Iyy, g, mu_g});
    
    % Constraint Jacobian
    
    % Generate constraint jacobian function
    fprintf('Constraint Jacobian function Gen...\n');
    matlabFunction(constraint_jacobian_nz,...
        'file','Jump2D/Jump2DConstraintJacobian',...
        'vars',{states, inputs, c_states, c_states0, dt, r_foot, m, Iyy, mu_g});
    
    % Generate constraint jacobian sparsity function
    fprintf('Constraint Jacobian Sparsity function Gen...\n')
    matlabFunction(row_index_CJ, col_index_CJ,...
        'file','Jump2D/Jump2DConstraintJacobianSP',...
        'vars', {iter, NUM_X, NUM_C});
    
    % Generate final constraint jacobian function
    fprintf('Final Constraint Jacobian function Gen...\n');
    matlabFunction(constraint_jacobian_final_nz,...
        'file','Jump2D/Jump2DConstraintJacobianFinal',...
        'vars',{states, inputs, c_states, c_states0, dt, r_foot, m, Iyy, mu_g});
    
    % Generate final constraint jacobian sparsity function
    fprintf('Final Constraint Jacobian Sparsity function Gen...\n')
    matlabFunction(row_index_final_CJ, col_index_final_CJ,...
        'file','Jump2D/Jump2DConstraintJacobianFinalSP',...
        'vars', {iter, NUM_X, NUM_C});
    
    % Generate initial constraint jacobian function
    fprintf('Final Constraint Jacobian function Gen...\n');
    matlabFunction(constraint_jacobian_initial_nz,...
        'file','Jump2D/Jump2DConstraintJacobianInitial',...
        'vars',{states, inputs, c_states, c_states0, dt, r_foot, m, Iyy, mu_g});
    
    % Generate initial constraint jacobian sparsity function
    fprintf('Initial Constraint Jacobian Sparsity function Gen...\n')
    matlabFunction(row_index_initial_CJ, col_index_initial_CJ,...
        'file','Jump2D/Jump2DConstraintJacobianInitialSP',...
        'vars', {iter, NUM_X, NUM_C});
    
    % Hessian
    
    % Generate lagrangian hessian function
    fprintf('Lagrangian Hessian function Gen...\n')
    matlabFunction(lagrangian_hessian_nz,...
        'file','Jump2D/Jump2DLagrangianHessian',...
        'vars',{c_states, dt, Q, R, obj_factor, lambda, Iyy});
    
    % Generate hessian sparsity function
    fprintf('Lagrangian Hessian Sparsity function Gen...\n')
    matlabFunction(row_index_H, col_index_H,...
        'file','Jump2D/Jump2DLagrangianHessianSP',...
        'vars', {iter, NUM_X});
    
    % Generate final lagrangian hessian function
    fprintf('Final Lagrangian Hessian function Gen...\n')
    matlabFunction(lagrangian_hessian_final_nz,...
        'file','Jump2D/Jump2DLagrangianHessianFinal',...
        'vars',{c_states, dt, Q, R, obj_factor, lambda, Iyy});
    
    % Generate final hessian sparsity function
    fprintf('Final Lagrangian Hessian Sparsity function Gen...\n')
    matlabFunction(row_index_final_H, col_index_final_H,...
        'file','Jump2D/Jump2DLagrangianHessianFinalSP',...
        'vars', {iter, NUM_X});
    
    % Generate initial lagrangian hessian function
    fprintf('Initial Lagrangian Hessian function Gen...\n')
    matlabFunction(lagrangian_hessian_initial_nz,...
        'file','Jump2D/Jump2DLagrangianHessianInitial',...
        'vars',{c_states, dt, Q, R, obj_factor, lambda, Iyy});
    
    % Generate initial hessian sparsity function
    fprintf('Initial Lagrangian Hessian Sparsity function Gen...\n')
    matlabFunction(row_index_final_H, col_index_final_H,...
        'file','Jump2D/Jump2DLagrangianHessianInitialSP',...
        'vars', {iter, NUM_X});
    
    fprintf('DONE generating MATLAB functions: %f s\n\n',toc(functions_timer));
end

% Finalize generation statistics
fprintf(['DONE generating symbolic functions!\n',...
    '   Total time taken: %f s\n\n'],toc(SymbolicJump2D_timer))