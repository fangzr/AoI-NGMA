function [pre,op] = Get_data_L_EE(result,result_optimal)
len = length(result);
op = zeros(1,len);
pre = zeros(1,len);


%% EE_A_bitrate

for i=1:len
    pre(1,i)=10 * log10((result(i).E_con)*1e3);
    op(1,i)=10 * log10((result_optimal(i).E_con)*1e3) ;
end