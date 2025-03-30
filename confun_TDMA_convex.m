function [c,ceq] = confun_TDMA_convex(x,model)
% Nonlinear inequality constraints
T_p = model.T_p;
Num_MTCD = model.Num_MTCD;
L=model.L;
B=model.B;
gamma=model.gamma(Num_MTCD);

% c = [(x(1)-T_p) * ( 2 ^ (N * L/((x(1)-T_p)*B))-1)/gamma-x(4)*T_p;
% 	0-lambda_temp*x(1);
% 	lambda_temp*x(1)-1;
%     ( 2 ^ (N * L/((x(1)-T_p)*B))-1)/gamma-model.P_max
% 	];

c = [
    (x(1)-T_p)/Num_MTCD * ( 2 ^ (Num_MTCD * L/((x(1)-T_p)*B))-1)/gamma-x(3)*T_p; % 31a
    ( 2 ^ (Num_MTCD * L/((x(1)-T_p)*B))-1)/gamma-model.P_max;%31a,新加的功率约束 ( 2 ^ (Num_MTCD * L/((x(1)-T_p)*B))-1)/gamma-model.P_max;
     model.EE * (x(1)-T_p)* (2 ^ (Num_MTCD * L/((x(1)-T_p)*B))-1)- Num_MTCD * L * gamma %30
	];

% Nonlinear equality constraints
ceq = [];
end