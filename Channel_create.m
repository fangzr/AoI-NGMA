%参数初始化后，修改信道模型

function model = Channel_create(model,num)

Num_MTCD = model.Num_MTCD;
 D_max = 5;
 %最小节点距离 (m)
 D_min = 3;
  % 距离集合
 D_PM = linspace(D_min,D_max,Num_MTCD);

 if(num==1)
    %% TDMA

    % 信道特性

    B = 5e6;%信道带宽（Hz）
    h = 1e-3 * 1 * D_PM.^(-2);%EH到MTCD衰减系数
    g = 1e-3 * 1 * D_PM.^(-2);%MTCD到BS衰减系数
    N = 1e-3 * 10^(-11) * B;%噪声功率 -60dbm/Hz
    gamma = (g.^2)/N;




elseif num == 2
    %% FDMA


    % 信道特性
    
    B = 5e6/Num_MTCD;%信道带宽（Hz）
%     D_p = 3;%EH到MTCD距离（m）
%     D_s = 3;%MTCD到BS距离（m）

    h = 1e-3 * 1 * D_PM.^(-2);%EH到MTCD衰减系数
    g = 1e-3 * 1 * D_PM.^(-2);%MTCD到BS衰减系数
    N = 1e-3 * 10^(-11) * B;%噪声功率 -60dbm/Hz
    gamma = (g.^2)/N;


    
else
    %% NOMA


    
    B = 5e6;%信道带宽（Hz）


    h = 1e-3 * 1 * D_PM.^(-2);%EH到MTCD衰减系数
    g = 1e-3 * 1 * D_PM.^(-2);%MTCD到BS衰减系数
    N = 1e-3 * 10^(-11) * B;%噪声功率 -60dbm/Hz
    gamma = (g.^2)/N;
    
    
    
end



model.h=h;
model.g=g;
model.gamma=gamma;

