clc
clear

%% Initialization
model = CreateModel(2);
L_vector = model.L_vector;
Num_iter = size(L_vector,2);
% Create Empty Particle Structure
empty_result.Lmabda=[];
empty_result.T_b=[];
empty_result.M=[];
empty_result.P_r=[];
empty_result.AoI=[];
empty_result.EE=[];
empty_result.P_tr=[];
empty_result.E_con=[];
empty_result.AoI_no_vac=[];
empty_result.E_con_no_vac=[];
% Create Result Matrix
result=repmat(empty_result,Num_iter,1);

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



for k=1:Num_iter
    x0 = [1,0.012,4,1.57];%initialize the first point
%     optimal_x0 = [10,0.0165791417420164,4.00000799994102,1.56505191158835;11,0.0165791339902409,4.00000879989178,1.56539671596784;12,0.0165791261364347,4.00000959981151,1.56552527500388;13,0.0165791181762074,4.00001039967906,1.56557153348247;14,0.0165791102296303,4.00001119999224,1.56580983008450;15,0.0165791020546754,4.00001199986563,1.57002380645611;16,0.0165790939324919,4.00001279981783,1.56674517557702;17,0.0165790856754498,4.00001359988643,1.56680808979727;18,0.0165790773323568,4.00001439992155,1.56815200467740;19,0.0165790689747012,4.00001519998567,1.56582320115325;20,0.0165790605209443,4.00001599995101,1.56648797728109];
%     x0=optimal_x0(k,:);
    last_x=x0;
    now_x=x0;
    i=1;
    model.L = L_vector(k);
    result_AoI = zeros(1,Num_iter+1);
    temp_norm = zeros(1,Num_iter+1);
    % temp_norm(i)=norm(last_x-x);
    AoI_optimal = 1e2;
    %% CCP 
    while (i<10)
    % Input the last optimal solution x_(k)
    i = i+1;
    last_x=now_x;
    eta1 = (last_x(2))^2/(2 * (1-last_x(1) * last_x(2))^2) - last_x(3)/(2 * (last_x(1))^2);

    eta2 = 1/(2 * (1-last_x(1) * last_x(2))^2)-0.5;

    eta3 = 1/(2 * last_x(1));

    fun = @(x) x(2) + 1/(2 * x(1)) + x(1) * eta1 + x(2) * eta2 + x(3) * eta3;
    % %3变量
    % A = [-1,0,0;
    % 	  1,0,0;
    % 	  0,0,1;
    % 	  0,0,-1;
    % 	  0,-1,0;
    % 	  0,1,0
    % 	];
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
          model.M_max;
          -1 * model.M_min;
          -1*model.T_b_min;
          model.T_b_max;
          model.E_b_max/model.T_p;
          1+last_x(1)*last_x(2) 
        ];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    c_confun = @(x)confun_FDMA(x,model);
    [now_x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(x)confun(x,model));
    result_AoI_temp = now_x(2) + now_x(1)* ((now_x(2))^2) /(2 * (1 - now_x(1) * now_x(2)))  +( now_x(3)+1)/(2 * now_x(1));
    if exitflag>0  && result_AoI_temp>0 && result_AoI_temp<AoI_optimal
        AoI_optimal = result_AoI_temp;
        temp_norm(i)=norm(last_x-now_x);
        % Obtain the final result in the i-th iteration
        result(k).Lmabda = now_x(1);
        result(k).T_b = now_x(2);
        result(k).M = now_x(3);
        result(k).P_r = now_x(4);
        result(k).AoI =  result_AoI_temp;
        result(k).EE = model.gamma(Num_MTCD) * model.Num_MTCD *  model.L / ((now_x(2)-model.T_p) * (2 ^ (model.Num_MTCD *  model.L/(model.B * (now_x(2)-model.T_p)))-1));
        result(k).P_tr = (2 ^ ( model.L/((now_x(2)-model.T_p)*model.B))-1)/model.gamma(Num_MTCD);
        result(k).E_con = (result(k).P_tr *(result(k).T_b-model.T_p)/(result(k).T_b)+model.P_work) * result(k).Lmabda * result(k).T_b + (1-result(k).Lmabda * result(k).T_b) * model.P_sleep;
         %无休眠结果
        result(k).AoI_no_vac =  0.5*(now_x(2)) + 1/now_x(1) + now_x(2)/(2*(1-now_x(1)*now_x(2)));
        result(k).E_con_no_vac =(result(k).P_tr * (result(k).T_b-model.T_p)/(result(k).T_b) +model.P_work) * result(k).Lmabda * result(k).T_b + (1-result(k).Lmabda * result(k).T_b) * model.P_work;
    end
    end


 %% Find the average value of energy consumption of nodes at different distances
	for u=1:MTCD_num
		
		P_tr_average(k,u)= (2 ^ (model.L/((result(k).T_b-model.T_p)*model.B))-1)/model.gamma(u);%注意和TDMA的不同！！！
		
		E_con_average(k,u) = (P_tr_average(k,u) *(result(k).T_b-model.T_p)/(result(k).T_b )+model.P_work) * result(k).Lmabda * result(k).T_b + (1-result(k).Lmabda * result(k).T_b) *  model.P_sleep;
		
 		E_con_no_vac_average(k,u) =(P_tr_average(k,u) * (result(k).T_b-model.T_p)/(result(k).T_b ) +model.P_work) * result(k).Lmabda * result(k).T_b + (1-result(k).Lmabda * result(k).T_b) * model.P_work;
	end

end



P_tr_temp = mean(P_tr_average,2);
E_con_temp = mean(E_con_average,2);
E_con_no_vac_temp =  mean(E_con_no_vac_average,2);
    
for k=1:Num_iter
    result_temp(k).P_tr = P_tr_temp(k);
    result_temp(k).E_con = E_con_temp(k);
    result_temp(k).E_con_no_vac = E_con_no_vac_temp(k);
end

%% FDMA Algorithm 2 for convex optimization

empty_result.Lmabda=[];
empty_result.T_b=[];
empty_result.M=[];
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
         x0 = [0.0215099768994962,2.00000064000003,2.02319463472494];%initialize the first point
    fun = @(x) x(1) + lambda_temp * (x(1))^2/( 2 * (1 - lambda_temp * x(1))) + (x(2)+1)/(2 * lambda_temp);


    A = [ 
          0,1,0;
          0,-1,0;
          -1,0,0;
          1,0,0;
          0,0,1;
          lambda_temp,0,0
        ];

    b = [
          model.M_max;
          -1 * model.M_min;
          -1*model.T_b_min;
          model.T_b_max;%31
          model.E_b_max/model.T_p;
          1
        ];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    c_confun = @(x)confun_FDMA_convex(x,model);
    [x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(x)confun_FDMA_convex(x,model));
     AoI_temp =  x(1) + lambda_temp * (x(1))^2/( 2 * (1 - lambda_temp * x(1))) + (x(2)+1)/(2 * lambda_temp);
    if AoI_temp<local_optimal_lambda && exitflag>0
        result_optimal(k).Lmabda = lambda_temp;
        result_optimal(k).T_b = x(1);
        result_optimal(k).M = x(2);
        result_optimal(k).P_r = x(3);
        result_optimal(k).AoI =  AoI_temp;
        result_optimal(k).P_tr = (2 ^ (model.L/((result_optimal(k).T_b-model.T_p)*model.B))-1)/model.gamma(Num_MTCD);
        result_optimal(k).E_con = (result_optimal(k).P_tr * (result_optimal(k).T_b-model.T_p)/(result_optimal(k).T_b) +model.P_work) * result_optimal(k).Lmabda * result_optimal(k).T_b + (1-result_optimal(k).Lmabda * result_optimal(k).T_b) * model.P_sleep;%Note that TDMA consumes transmission power only in the corresponding time slot
        local_optimal_lambda = AoI_temp;
    end
    end
	%% Find the average value of energy consumption of nodes at different distances
	for u=1:MTCD_num
		
		P_tr_average_optimal(k,u)= (2 ^ ( model.L/((result_optimal(k).T_b-model.T_p)*model.B))-1)/model.gamma(u);
		
		E_con_average_optimal(k,u) = (P_tr_average_optimal(k,u) * (result_optimal(k).T_b-model.T_p)/(result_optimal(k).T_b) +model.P_work) * result_optimal(k).Lmabda * result_optimal(k).T_b + (1-result_optimal(k).Lmabda * result_optimal(k).T_b) *  model.P_sleep;
		
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
save('E:\IEEE_LaTex\21-JSAC\Simulation\Data\FDMA_ST_L.mat');