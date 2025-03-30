%% Simulation Parameter Settings

function model = CreateModel(num)

%% Simulation parameters for Algorithm 2
model.step_num = 20;
model.L_vector = [60:20:160];
model.Lambda_max = [10:2:20];
model.Num_MTCDseq = [10:1:16];
gamma = zeros(1,100);
P_max_cons = 0.4; % Transmit power constraint for MTCD
E_b_max_cons = 4e-2; % Received power at the energy station for energy transmission
M_max = 20;
M_min = 2;
CCP_lambda_min = 1;
CCP_lambda_max = 15;
T_b_max = 5e-2;
% Number of nodes
Num_MTCD = 10;
% Maximum node distance (m)
D_max = 5;
% Minimum node distance (m)
D_min = 3;
% Distance set
D_PM = linspace(D_min, D_max, Num_MTCD);
% Initialization
D_p = D_min;
D_s = D_p;

% %% Temporary setting for Figure 5, originally default was 80
L = 100;
% Channel gap
channel_J = 0.05;

%% Power settings

E_b_max = E_b_max_cons; % Charging energy (J)

P_max = P_max_cons; % Maximum transmit power of MTCD

P_work = 10e-2; % Energy consumption in active (working) state

P_sleep = 1e-2; % Energy consumption in sleep state

P_sc = 10e-2; % Switching cost energy consumption

EE = 1.4e5; % Energy efficiency
 
if num == 1

    %% TDMA

    % Channel characteristics

    B = 5e6; % Channel bandwidth (Hz)
%     D_p = 3; % Distance from EH to MTCD (m)
%     D_s = 3; % Distance from MTCD to BS (m)

    h = 1e-3 * 1 * D_PM.^(-2); % Path loss coefficient from EH to MTCD
    g = 1e-3 * 1 * D_PM.^(-2); % Path loss coefficient from MTCD to BS
    N = 1e-3 * 10^(-11) * B; % Noise power (-60 dBm/Hz)
    gamma = (g.^2) / N;

    EH_EE = 0.9; % Charging efficiency
    % Time slot constraints
    T_p = 1e-2; % Length of one frame transmission/charging slot (s)
    T_b_min = T_p + 5e-3;
    T_s_max = 0.5;
    T_s_min = 1e-2;
    
    % Initial value
%     x0_MV = 
elseif num == 2
    %% FDMA

    % Channel characteristics
    
    B = 5e6 / Num_MTCD; % Channel bandwidth (Hz)
%     D_p = 3; % Distance from EH to MTCD (m)
%     D_s = 3; % Distance from MTCD to BS (m)

    h = 1e-3 * 1 * D_PM.^(-2); % Path loss coefficient from EH to MTCD
    g = 1e-3 * 1 * D_PM.^(-2); % Path loss coefficient from MTCD to BS
    N = 1e-3 * 10^(-11) * B; % Noise power (-60 dBm/Hz)
    gamma = (g.^2) / N;

    EH_EE = 0.9; % Charging efficiency
    % Time slot constraints
    T_p = 1e-2; % Length of one frame transmission/charging slot (s)
    T_b_min = T_p + 5e-3;
    T_s_max = 0.5;
    T_s_min = 10e-3;

else
    %% NOMA

    B = 5e6; % Channel bandwidth (Hz)
%     D_p = 3; % Distance from EH to MTCD (m)
%     D_s = 3; % Distance from MTCD to BS (m)

    h = 1e-3 * 1 * D_PM.^(-2); % Path loss coefficient from EH to MTCD
    g = 1e-3 * 1 * D_PM.^(-2); % Path loss coefficient from MTCD to BS
    N = 1e-3 * 10^(-11) * B; % Noise power (-60 dBm/Hz)
    gamma = (g.^2) / N;

    EH_EE = 0.9; % Charging efficiency
    % Time slot constraints
    T_p = 1e-2; % Length of one frame transmission/charging slot (s)
    T_b_min = T_p + 5e-3;
    T_s_max = 0.5;
    T_s_min = 10e-3;

    % Power constraints
    E_b_max = E_b_max_cons * Num_MTCD; % Charging energy (J)
    P_max = P_max_cons * Num_MTCD; % Maximum transmit power of MTCD
    P_work = 1e-1; % Energy consumption in active (working) state
    P_sleep = 1e-2; % Energy consumption in sleep state
    EE = 1.4e5; % Energy efficiency
end

model.B = B;
% model.D_p = D_p;
% model.D_s = D_s;
model.h = h;
model.g = g;
model.N = N;
model.gamma = gamma;
model.Num_MTCD = Num_MTCD;
model.L = L;
model.EH_EE = EH_EE;
model.T_p = T_p;
model.T_b_min = T_b_min;
model.T_b_max = T_b_max;
model.CCP_lambda_min = CCP_lambda_min;
model.CCP_lambda_max = CCP_lambda_max;
model.E_b_max = E_b_max;
model.T_s_max = T_s_max;
model.P_max = P_max;
model.P_work = P_work;
model.P_sc = P_sc;
model.P_sleep = P_sleep;
model.EE = EE;
model.T_s_min = T_s_min;
model.M_max = M_max;
model.M_min = M_min;
model.channel_J = channel_J; % Channel gap

% Charging rate threshold
model.Recharge_rate_threshold = E_b_max / T_p;

end
