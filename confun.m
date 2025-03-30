function [c,ceq] = confun(x,model)
% Nonlinear inequality constraints
T_p = model.T_p;
Num_MTCD = model.Num_MTCD;
L=model.L;
B=model.B;
gamma=model.gamma(Num_MTCD);
%% Consider EH power
% c = [(x(2)-T_p) * ( 2 ^ (N * L/((x(2)-T_p)*B))-1)/gamma-x(4)*T_p;
% 	0-x(1)*x(2);
% 	x(1)*x(2)-1;
%     ( 2 ^ (N * L/((x(2)-T_p)*B))-1)/gamma-model.P_max
% 	];
%% Don't Consider EH power
c = [
    (x(2)-T_p)/Num_MTCD * ( 2 ^ (Num_MTCD * L/((x(2)-T_p)*B))-1)/gamma-x(4)*T_p;
    0-x(1)*x(2);
    ( 2 ^ (Num_MTCD * L/((x(2)-T_p)*B))-1)/gamma-model.P_max;
     model.EE * (x(2)-T_p)* (2 ^ (Num_MTCD * L/((x(2)-T_p)*B))-1)- Num_MTCD * L * gamma
	];

% Nonlinear equality constraints
ceq = [];
end