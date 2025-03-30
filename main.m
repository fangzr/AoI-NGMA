clc
clear
close 
addpath Data figures Main-MV-ST-NVaca Main-MV-ST Main-MV;
%% Run the main program Figure 1-4
% Main_TDMA_ST
% clear
% Main_TDMA_MV
% clear
% Main_FDMA_ST
% clear
% Main_FDMA_MV
% clear
% Main_NOMA_ST
% clear
% Main_NOMA_MV
% 
%% Draw
% figure(1)
% Plot_fig_AoI_L_ST(CreateModel(1));
% figure(2)
% Plot_fig_AoI_L(CreateModel(1));
% figure(3)
% Plot_fig_EE_L_ST(CreateModel(1));
% figure(4)
% Plot_fig_EE_L_MV(CreateModel(1));
%% Figure 5-6
% 
% Main_TDMA_ST_lambda
% clear
% Main_TDMA_MV_lambda
% clear
% Main_FDMA_ST_lambda
% clear
% Main_FDMA_MV_lambda
% clear
% Main_NOMA_ST_lambda
% clear
% Main_NOMA_MV_lambda
% 
% Plot_fig_AoI_Lambda(CreateModel(1));


%% PP AoI vs L
% clear
% figure(5)
% Plot_fig_pp_AoI_L(CreateModel(1));

%% Figure 7


Main_TDMA_MV_Num
clear

Main_FDMA_MV_Num
clear

Main_NOMA_MV_Num
Plot_fig_EE_AoI_Num(CreateModel(1));


%% finish
f = warndlg('Simulation Done！','Warning');