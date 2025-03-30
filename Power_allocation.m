%% Calculate each power allocation in NOMA scheme
function [P_NOMA,P_average,E_con] = Power_allocation(result,model)
Num = model.Num_MTCD;
channel_J = model.channel_J;
P_NOMA_temp = zeros(1,Num);
Cap = (channel_J+1) * model.L/(result.T_b - model.T_p);

temp = 2^(Cap/model.B)-1;

%Turn the gain vector, become the first largest, the Nth smallest, first solve the Nth SINR
gamma_rev =  fliplr(model.gamma);

for i=1:Num
    j=Num-i+1;
    if i==1
         P_NOMA_temp(j)=temp/gamma_rev(i);
    else
        P_NOMA_temp(j)= temp * ( 1+ sum(gamma_rev(j+1:Num) .*(P_NOMA_temp(j+1:Num))))/gamma_rev(j);
    end
end
%Turn it back again to get the node power distribution, the furthest one gets the most power should be the biggest
 P_NOMA = fliplr(P_NOMA_temp);
 
 P_average = sum(P_NOMA,2)/Num;
 
 E_con = (P_average * (result.T_b-model.T_p)/(result.T_b) +model.P_work) * result.Lmabda * result.T_b + (1-result.Lmabda * result.T_b)* (0.9 * model.P_sleep + 0.1 * model.P_sc );