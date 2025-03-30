clc
clear

%% Initialization
model = CreateModel(1);
L_vector = model.L_vector;
Num_iter = size(L_vector,2);
% Create Empty Particle Structure
empty_result.Lmabda=[];
empty_result.T_b=[];
empty_result.T_s=[];
empty_result.P_r=[];
empty_result.AoI=[];
empty_result.EE=[];
empty_result.P_tr=[];
empty_result.E_con=[];
empty_result.AoI_no_vac=[];
empty_result.E_con_no_vac=[];

% Create Result Matrix
result=repmat(empty_result,Num_iter,1);

%% Add Var
result_temp=repmat(empty_result,Num_iter,1);
EE_average = zeros(Num_iter,model.Num_MTCD);
P_tr_average = zeros(Num_iter,model.Num_MTCD);
E_con_average = zeros(Num_iter,model.Num_MTCD);
E_con_no_vac_average =zeros(Num_iter,model.Num_MTCD);
P_tr_average_optimal = zeros(Num_iter,model.Num_MTCD);
E_con_average_optimal = zeros(Num_iter,model.Num_MTCD);
MTCD_num = model.Num_MTCD;
Num_MTCD = model.Num_MTCD;
result_temp_optimal=repmat(empty_result,Num_iter,1);
%% TDMA-CCP

for k=1:Num_iter
    %%initialize the first point
    %(packet arrival rate,  busy time slot, sleep-scheduling slot, received energy power)
    x0 = [0.05,0.15,0.4,2.5];
    last_x=x0;
    now_x=x0;
    i=1;
    model.L = L_vector(k); %修改model
    result_AoI = zeros(1,Num_iter+1);
    temp_norm = zeros(1,Num_iter+1);
    % temp_norm(i)=norm(last_x-x);
    j=1;
    AoI_optimal = 1e2;

    %% CCP 
    while (i<20)
    % Input the last optimal solution x_(k)
    i = i+1;
    last_x=now_x;
    eta1 = (last_x(2))^2 / (2 * (1-last_x(1) * last_x(2))^2);
    eta2 = (last_x(2))^2 /2;

    fun = @(x) (x(3)+x(2))/2 + 1/(x(1)) + eta1 * ( eta2 * x(1) + 0.5 * x(2));

    A = [-1,0,0,0;
          1,0,0,0;
          0,0,1,0;
          0,0,-1,0;
          0,-1,0,0;
          0,1,0,0;
          0,0,0,1;
          last_x(2),last_x(1),0,0
        ];

    b = [-1*model.CCP_lambda_min;
          model.CCP_lambda_max;
          model.T_s_max;
          -1 * model.T_s_min;
          -1*model.T_b_min;
          model.T_b_max;
          model.E_b_max/model.T_p;
          1+last_x(1)*last_x(2)
        ];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    c_confun = @(x)confun(x,model);
%     option = optimoptions('fmincon','Algorithm','sqp');
    [now_x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(x)confun(x,model));
    result_AoI_temp = 0.5*(now_x(2)+now_x(3)) + 1/now_x(1) + now_x(2)/(2*(1-now_x(1)*now_x(2)));
        if exitflag>0 && result_AoI_temp<AoI_optimal && result_AoI_temp>0
            AoI_optimal = result_AoI_temp;
            temp_norm(i)=norm(last_x-now_x);

            result_AoI(j) = result_AoI_temp;j=j+1;
            % Obtain the final result in the i-th iteration
            result(k).Lmabda = now_x(1);
            result(k).T_b = now_x(2);
            result(k).T_s = now_x(3);
            result(k).P_r = now_x(4);
            result(k).AoI =  result_AoI_temp;
			
            result(k).EE = model.gamma(MTCD_num) * model.Num_MTCD *  model.L / ((now_x(2)-model.T_p) * (2 ^ ((model.channel_J + 1) * model.Num_MTCD *  model.L/(model.B * (now_x(2)-model.T_p)))-1));
			
            result(k).P_tr = (2 ^ ((model.channel_J + 1) * model.Num_MTCD * model.L/((now_x(2)-model.T_p)*model.B))-1)/model.gamma(MTCD_num);
			
            result(k).E_con = (result(k).P_tr * (result(k).T_b-model.T_p)/(result(k).T_b * model.Num_MTCD) +model.P_work) * result(k).Lmabda * result(k).T_b + (1-result(k).Lmabda * result(k).T_b) * model.P_sleep;%Note that TDMA consumes transmission power only in the corresponding time slot
            %No sleep
            result(k).AoI_no_vac =  0.5*(now_x(2)) + 1/now_x(1) + now_x(2)/(2*(1-now_x(1)*now_x(2)));
            result(k).E_con_no_vac =(result(k).P_tr * (result(k).T_b-model.T_p)/(result(k).T_b * model.Num_MTCD) +model.P_work) * result(k).Lmabda * result(k).T_b + (1-result(k).Lmabda * result(k).T_b) * model.P_work;
        end
    end
    %% Find the average value of energy consumption of nodes at different distances
	for u=1:MTCD_num
		EE_average(k,u)= model.gamma(u) * model.Num_MTCD *  model.L / ((result(k).T_b-model.T_p) * (2 ^ ((model.channel_J + 1) * model.Num_MTCD *  model.L/(model.B * (result(k).T_b-model.T_p)))-1));
		
		P_tr_average(k,u)= (2 ^ ((model.channel_J + 1) * model.Num_MTCD * model.L/((result(k).T_b-model.T_p)*model.B))-1)/model.gamma(u);
		
		E_con_average(k,u) = (P_tr_average(k,u) * (result(k).T_b-model.T_p)/(result(k).T_b ) +model.P_work) * result(k).Lmabda * result(k).T_b + (1-result(k).Lmabda * result(k).T_b) * (0.9 * model.P_sleep + 0.1 * model.P_sc );
		
		E_con_no_vac_average(k,u) =(P_tr_average(k,u) * (result(k).T_b-model.T_p)/(result(k).T_b) +model.P_work) * result(k).Lmabda * result(k).T_b + (1-result(k).Lmabda * result(k).T_b) * model.P_work;
	end
end

EE_temp =mean(EE_average,2);
P_tr_temp = mean(P_tr_average,2);
E_con_temp = mean(E_con_average,2);
E_con_no_vac_temp =  mean(E_con_no_vac_average,2);
    
for k=1:Num_iter
    result_temp(k).EE = EE_temp(k);
    result_temp(k).P_tr = P_tr_temp(k);
    result_temp(k).E_con = E_con_temp(k);
    result_temp(k).E_con_no_vac = E_con_no_vac_temp(k);
end

%% TDMA Algorithm 2 for convex optimization
% Create Empty Particle Structure
empty_result.Lmabda=[];
empty_result.T_b=[];
empty_result.T_s=[];
empty_result.P_r=[];
empty_result.AoI=[];
empty_result.EE=[];
empty_result.P_tr=[];
empty_result.E_con=[];

% Create Result Matrix
result_optimal=repmat(empty_result,Num_iter,1);

for k=1:Num_iter
    model.L = L_vector(k);
    Lambda_seq = [model.CCP_lambda_min:(model.CCP_lambda_max-model.CCP_lambda_min)/(model.step_num-1):model.CCP_lambda_max];
    local_optimal_lambda = 1e2;
    
    
    for i=1:model.step_num
        lambda_temp = Lambda_seq(i);
         x0 = [0.03,0.02,0.5];%initialize the first point
    fun = @(x) (x(1)+x(2))/2 + 1/(lambda_temp) + x(1)/(2 * (1 - lambda_temp * x(1)));

    A = [ 
          0,1,0;
          0,-1,0;
          -1,0,0;
          1,0,0;
          0,0,1;
          lambda_temp,0,0
        ];

    b = [
          model.T_s_max;
          -1 * model.T_s_min;
          -1*model.T_b_min;
          model.T_b_max;%31
          model.E_b_max/model.T_p;
          1
        ];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    c_confun = @(x)confun_TDMA_convex(x,model);
    [x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(x)confun_TDMA_convex(x,model));
    AoI_temp =  0.5*(x(1)+x(2)) + 1/lambda_temp + x(1)/(2*(1-lambda_temp*x(1)));
    if AoI_temp<local_optimal_lambda && exitflag>0
        result_optimal(k).Lmabda = lambda_temp;
        result_optimal(k).T_b = x(1);
        result_optimal(k).T_s = x(2);
        result_optimal(k).P_r = x(3);
        result_optimal(k).AoI =  AoI_temp;
        result_optimal(k).P_tr = (2 ^ ((model.channel_J + 1) * model.Num_MTCD * model.L/((result_optimal(k).T_b-model.T_p)*model.B))-1)/model.gamma(MTCD_num);
        result_optimal(k).E_con = (result_optimal(k).P_tr * (result_optimal(k).T_b-model.T_p)/(result_optimal(k).T_b ) +model.P_work) *result_optimal(k).Lmabda * result_optimal(k).T_b + (1-result_optimal(k).Lmabda * result_optimal(k).T_b) * (0.9 * model.P_sleep + 0.1 * model.P_sc );%Note that TDMA consumes transmission power only in the corresponding time slot
        local_optimal_lambda = AoI_temp;
    end
    end
	%% Find the average value of energy consumption of nodes at different distances
	for u=1:MTCD_num
		
		P_tr_average_optimal(k,u)= (2 ^ ((model.channel_J + 1) * model.Num_MTCD * model.L/((result_optimal(k).T_b-model.T_p)*model.B))-1)/model.gamma(u);
		
		E_con_average_optimal(k,u) = (P_tr_average_optimal(k,u) * (result_optimal(k).T_b-model.T_p)/(result_optimal(k).T_b ) +model.P_work) * result_optimal(k).Lmabda * result_optimal(k).T_b + (1-result_optimal(k).Lmabda * result_optimal(k).T_b) * (0.9 * model.P_sleep + 0.1 * model.P_sc );
		
	end
end

% EE_temp =mean(EE_average,2);
P_tr_temp = mean(P_tr_average_optimal,2);
E_con_temp = mean(E_con_average_optimal,2);

    
for k=1:Num_iter
%     result_temp_optimal(k).EE = EE_temp(k);
    result_temp_optimal(k).P_tr = P_tr_temp(k);
    result_temp_optimal(k).E_con = E_con_temp(k);
%     result_temp_optimal(k).E_con_no_vac = E_con_no_vac_temp(k);
end
save('E:\IEEE_LaTex\21-JSAC\Simulation\Data\TDMA_MV_L.mat');
