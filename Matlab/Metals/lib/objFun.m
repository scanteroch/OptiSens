% Function: Entropy-based convex objective function
function [f, g] = objFun(candSol,NmaxSns,NmaxAct,N_smpl,Nvar,x_cost,...
    y_cost,P_sns_eval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input paramters:
%       candSol         -   candidate solution
%       NmaxAct         -   Maximum number of actuators
%       NmaxSns         -   Maximum number of sensors
%       N_smpl          -   Number of model parameter samples to evaluate
%                           the objective function
%       N_var           -   Number of model parameters
%       x_cost, y_cost  -   X and Y coordinates of the interpolating points
%                           used to build the cost function (by pchip())
%       P_sns_eval      -   Model-gradient multiplication matrices
%
%   Output parameters:
%       f               -   Objective function evaluation
%       g               -   Gradient eval. of the the objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z=candSol(1:NmaxAct);
    w=candSol(NmaxAct+1:end-1);
    n_sns_act=candSol(end);
    
    f_aux=zeros(1,N_smpl);
    g_aux1=zeros(NmaxAct,N_smpl);
    g_aux2=zeros(NmaxSns,N_smpl);
    for i=1:N_smpl
        P_aux=P_sns_eval(:,:,:,:,i);
        Q=permute(sum(reshape(z.*reshape(permute (sum(reshape(w.*reshape(permute (P_aux,[3 1 2 4]), NmaxAct,Nvar*Nvar*NmaxSns),NmaxAct,Nvar,Nvar,NmaxSns),1),[4,2,3,1]),NmaxSns,Nvar*Nvar),NmaxAct,Nvar,Nvar),1),[2,3,1]);
        f_aux(i)=log(det(Q));
        
        for l=1:NmaxAct
            g_aux1(l,i)=trace(Q\sum(permute(reshape(w.*reshape(permute(P_aux(:,:,:,l),[3 1 2]),NmaxSns,Nvar*Nvar),NmaxSns,Nvar,Nvar),[2,3,1]),3));
        end

        for j=1:NmaxSns
            g_aux2(j,i)=trace(Q\sum(permute(reshape(z.*reshape(permute(permute(P_aux(:,:,j,:),[1,2,4,3]),[3 1 2]),NmaxSns,Nvar*Nvar),NmaxSns,Nvar,Nvar),[2,3,1]),3));
        end

    end
    
    % Cost function evaluation
    g_eval=pchip(x_cost,y_cost,n_sns_act);
    
    % Monte Carlo approximation of the objective function:    
    f_MC = -1/N_smpl * sum(f_aux);
    f = f_MC + abs(f_MC) * g_eval;
    
    % Derivative of the cost function
    grdPchip=n_sns_act-1:0.01:n_sns_act+1;
    dgAux=pchip(x_cost,y_cost,grdPchip);
    slp= diff(dgAux) ./ diff(grdPchip);
    slpGrd=grdPchip(1:end-1);
    dgEval=interp1(slpGrd,slp,n_sns_act, 'PCHIP');
    
    % Jacobian (gradient) 
    g = zeros(NmaxSns+NmaxAct+1,1);
    g(1:end-1) = [-1/N_smpl.* sum(g_aux1,2)  -1/N_smpl.* sum(g_aux2,2)] * ...
        (1+g_eval.*f_MC./abs(f_MC));
    g(end) = abs(f_MC)*dgEval;
end
