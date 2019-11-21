
%% Evaluates the Linearized Dynamics
function Yk1 = DiscreteDynamics(Xk, Uk, System)

% Discrete State Space system
Xk1 = System.A*Xk + System.B*Uk + System.Dw*System.w;
Yk1 = System.C*Xk1;

end
