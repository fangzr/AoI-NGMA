clc
clear

%% Initialization
model = CreateModel(3);
L_vector = model.L_vector;
Num_iter = size(L_vector,2);
% Create Empty Particle Structure
empty_result.Lmabda=[];
empty_result.T_b=[];
empty_result.M=[];
empty_result.P_r=[];
empty_result.AoI=[];
% empty_result.EE=[];
empty_result.P_tr=zeros(1,model.Num_MTCD);
empty_result.E_con=[];
empty_result.P_tr_ave=[];
% Create Result Matrix
result=repmat(empty_result,Num_iter,1);



for k=1:Num_iter
    x0 = [10,0.03,7,2];%initialize the first point
    last_x=x0;
    x=x0;
    i=1;
    model.L = L_vector(k);
    result_AoI = zeros(1,Num_iter+1);
    temp_norm = zeros(1,Num_iter+1);
    % temp_norm(i)=norm(last_x-x);
    %% CCP 
    while (i<20)
    % Input the last optimal solution x_(k)
    i = i+1;
    last_x=x;
    eta1 = (last_x(2))^2/(2 * (1-last_x(1) * last_x(2))^2) - last_x(3)/(2 * (last_x(1))^2);

    eta2 = 1/(2 * (1-last_x(1) * last_x(2))^2)-0.5;

    eta3 = 1/(2 * last_x(1));

    fun = @(x) x(2) + 1/(2 * x(1)) + x(1) * eta1 + x(2) * eta2 + x(3) * eta3;
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
    c_confun = @(x)confun_NOMA(x,model);
    [x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(x)confun_NOMA(x,model));
        if exitflag>0
            result(k).Lmabda = x(1);
            result(k).T_b = x(2);
            result(k).M = x(3);
            result(k).P_r = x(4)/model.Num_MTCD;%Average EH energy per node
            result(k).AoI =  x(2) + 1/(2 * x(1)) + x(1)* ((x(2))^2) /(2 * (1 - x(1) * x(2)))  + x(3)/(2 * x(1));
            [result(k).P_tr,result(k).P_tr_ave,result(k).E_con] = Power_allocation(result(k),model);
            result(k).E_con = result(k).E_con -  (1-result(k).Lmabda * result(k).T_b)*((0.9 * model.P_sleep + 0.1 * model.P_sc )-model.P_sleep);
             %No sleep
            result(k).AoI_no_vac =  0.5*(x(2)) + 1/x(1) + x(2)/(2*(1-x(1)*x(2)));
            result(k).E_con_no_vac =(result(k).P_tr_ave * (result(k).T_b-model.T_p)/(result(k).T_b) +model.P_work) * result(k).Lmabda * result(k).T_b + (1-result(k).Lmabda * result(k).T_b) * model.P_work;
        end
    end

end


%% NOMA Algorithm 2 for convex optimization


% Create Result Matrix
result_optimal=repmat(empty_result,Num_iter,1);
for k=1:Num_iter
    model.L = L_vector(k);
    Lambda_seq = [model.CCP_lambda_min:(model.CCP_lambda_max-model.CCP_lambda_min)/(model.step_num-1):model.CCP_lambda_max];
    local_optimal_lambda = 1e2;
    for i=1:model.step_num
        lambda_temp = model.CCP_lambda_max;
         x0 = [0.15,2,5];%initialize the first point
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
    c_confun = @(x)confun_NOMA_convex(x,model);
    [x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(x)confun_NOMA_convex(x,model));
    AoI_temp =  x(1) + lambda_temp * (x(1))^2/( 2 * (1 - lambda_temp * x(1))) + (x(2)+1)/(2 * lambda_temp);
    if AoI_temp<=local_optimal_lambda && exitflag>0
        result_optimal(k).Lmabda = lambda_temp;
        result_optimal(k).T_b = x(1);
        result_optimal(k).M = x(2);
        result_optimal(k).P_r = x(3)/model.Num_MTCD;
        result_optimal(k).AoI =  AoI_temp;
        [result_optimal(k).P_tr,result_optimal(k).P_tr_ave,result_optimal(k).E_con] = Power_allocation(result_optimal(k),model);
        result_optimal(k).E_con = result_optimal(k).E_con -  (1-result_optimal(k).Lmabda * result_optimal(k).T_b)*((0.9 * model.P_sleep + 0.1 * model.P_sc )-model.P_sleep);
        local_optimal_lambda = AoI_temp;
    end
    end
end
save('E:\IEEE_LaTex\21-JSAC\Simulation\Data\NOMA_ST_L.mat');
