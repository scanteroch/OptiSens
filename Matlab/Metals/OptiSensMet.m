%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         OptiSens - Optimal Sensor and Actuator Placement
%                        ISOTROPIC PLATES
% Sergio Cantero Chinchilla
% V01 - 29/05/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 00 - Code initialisation 
restoredefaultpath
clearvars; close all; clc
addpath ./lib
softnm='OptiSens';
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'modal';
uiwait(msgbox(['\fontsize{9}',sprintf(['You will be guided through',...
    ' the optimisation software\n\n']),...
    '\it',sprintf('Note: enter the input variables as requested\n\n'),...
    '\rmClick OK to continue'],...
    [softnm,' - ISOTROPIC PLATES'],'help',CreateStruct));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 01 - Generate the input parameters for the optimisation code
[X_g_v,Y_g_v,X_d_v,Y_d_v,minDist,GridPtX,GridPtY,N_smpl,prior_smpl,...
    P_sns_eval,x_cost,y_cost,Nini]=inputPar(softnm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 02 - Sensor and actuator optimisation (./res/Opt_Sol.mat)
[x, X_g_v,Y_g_v,X_d_v,Y_d_v,minDist,GridPtX,GridPtY,prior_smpl,...
    NmaxSns,NmaxAct]=OptiSens(X_g_v,Y_g_v,X_d_v,Y_d_v,minDist,...
    GridPtX,GridPtY,N_smpl,prior_smpl,P_sns_eval,x_cost,y_cost,Nini);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 03 - Plot the resulting optimal layout (./res/Opt_Sol.pdf)
plotOASP(X_d_v,Y_d_v,X_g_v,Y_g_v,prior_smpl,GridPtX,GridPtY,...
    minDist,NmaxSns,NmaxAct,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 04 - Write the results in a text file (./res/Opt_Sol.txt)
writeOASP(x,NmaxAct,NmaxSns,GridPtX,GridPtY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Flag when finished
uiwait(msgbox(['\fontsize{9}',sprintf(['The optimisation of sensors and'...
    ' actuatros has finished.\n\n']),...
    '\it',sprintf('Note: The output files can be found in ./res\n\n'),...
    '\rmClick OK to continue'],...
    [softnm,' - ISOTROPIC PLATES'],'help',CreateStruct));