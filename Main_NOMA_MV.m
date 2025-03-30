clc
clear

%% Initialization
model = CreateModel(3);
L_vector = model.L_vector;
Num_iter = size(L_vector,2);
% Create Empty Particle Structure
empty_result.Lmabda=[];
empty_result.T_b=[];
empty_result.T_s=[];
empty_result.P_r=[];
empty_result.AoI=[];
% empty_result.EE=[];
empty_result.P_tr=zeros(1,model.Num_MTCD);
empty_result.E_con=[];
empty_result.P_tr_ave=[];
% Create Result Matrix
result=repmat(empty_result,Num_iter,1);
empty_result.AoI_no_vac=[];
empty_result.E_con_no_vac=[];

result_temp=repmat(empty_result,Num_iter,1);



for k=1:Num_iter
   x0 = [6.50000000000000,0.15,0.4,2.5];
    last_x=x0;
    x=x0;
    i=1;
    model.L = L_vector(k);
    result_AoI = zeros(1,Num_iter+1);
    temp_norm = zeros(1,Num_iter+1);
    % temp_norm(i)=norm(last_x-x);
    %% CCP 
    while (i<50)
    % Input the last optimal solution x_(k)
    i = i+1;
    last_x=x;
    eta1 = 1/((1-last_x(1)*last_x(2))^2);
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
    c_confun = @(x)confun_NOMA(x,model);
    [x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(x)confun(x,model));
    if exitflag>0
        temp_norm(i)=norm(last_x-x);
        result_AoI = 0.5*(x(2)+x(3)) + 1/x(1) + x(2)/(2*(1-x(1)*x(2)));
         % Obtain the final result in the i-th iteration
        result(k).Lmabda = x(1);
        result(k).T_b = x(2);
        result(k).T_s = x(3);
        result(k).P_r = x(4);
        result(k).AoI =  0.5*(x(2)+x(3)) + 1/x(1) + x(2)/(2*(1-x(1)*x(2)));
        % result(k).EE = model.gamma * model.Num_MTCD *  model.L / ((x(2)-model.T_p) * (2 ^ (model.Num_MTCD *  model.L/(model.B * (x(2)-model.T_p)))-1));
        [result(k).P_tr,result(k).P_tr_ave,result(k).E_con] = Power_allocation(result(k),model);
         %No sleep results
        result(k).AoI_no_vac =  0.5*(x(2)) + 1/x(1) + x(2)/(2*(1-x(1)*x(2)));
%         result(k).E_con_no_vac =(result(k).P_tr * (result(k).T_b-model.T_p)/(result(k).T_b * model.Num_MTCD) +model.P_work) * result(k).Lmabda * result(k).T_b +
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
        lambda_temp = Lambda_seq(i);
         x0 = [0.15,0.4,2.5];%initialize the first point
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
    c_confun = @(x)confun_NOMA_convex(x,model);
    [x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(x)confun_NOMA_convex(x,model));
    AoI_temp =  0.5*(x(1)+x(2)) + 1/lambda_temp + x(1)/(2*(1-lambda_temp*x(1)));
    if AoI_temp<local_optimal_lambda && exitflag>0
        result_optimal(k).Lmabda = lambda_temp;
        result_optimal(k).T_b = x(1);
        result_optimal(k).T_s = x(2);
        result_optimal(k).P_r = x(3)/model.Num_MTCD;
        result_optimal(k).AoI =  AoI_temp;
        [result_optimal(k).P_tr,result_optimal(k).P_tr_ave,result_optimal(k).E_con] = Power_allocation(result_optimal(k),model);
        local_optimal_lambda = AoI_temp;
    end
    end
end
save('E:\IEEE_LaTex\21-JSAC\Simulation\Data\NOMA_MV_L.mat');