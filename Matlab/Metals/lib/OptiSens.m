function [x, X_g_v,Y_g_v,X_d_v,Y_d_v,minDist,GridPtX,GridPtY,prior_smpl,...
    NmaxSns,NmaxAct]=OptiSens(X_g_v,Y_g_v,X_d_v,Y_d_v,minDist,GridPtX,...
    GridPtY,N_smpl,prior_smpl,P_sns_eval,x_cost,y_cost,Nini)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function that minimizes the convex entropy-based objective function.
%
%   Input paramters:
%       X_g_v, Y_g_v    -   X and Y coordinates of the vertices of the
%                           outer geometry of the plate
%       X_d_v, Y_d_v    -   X and Y coordinates of the vertices of the
%                           area of possible damage occurrence
%       minDist         -   Approximate distance between consecutive grid
%                           points (possible sensor/actuator locations)
%       GridPtX,GridPtY -   X and Y coordinates of the grid points of
%                           possible sensor/actuator locations
%       N_smpl          -   Number of model parameter samples to evaluate
%                           the objective function
%       prior_smpl      -   Set of N_smpl samples obtained from the prior
%                           distribution of the model parameters
%       P_sns_eval      -   Model-gradient multiplication matrices
%       x_cost, y_cost  -   X and Y coordinates of the interpolating points
%                           used to build the cost function (by pchip())
%       Nini            -   Guess solution for the optimal number of
%                           sensors
%
%   Output parameters:
%       x               -   Optimal solution containing z, w, and n
%       X_g_v, Y_g_v    -   X and Y coordinates of the vertices of the
%                           outer geometry of the plate
%       X_d_v, Y_d_v    -   X and Y coordinates of the vertices of the
%                           area of possible damage occurrence
%       minDist         -   Approximate distance between consecutive grid
%                           points (possible sensor/actuator locations)
%       GridPtX,GridPtY -   X and Y coordinates of the grid points of
%                           possible sensor/actuator locations
%       N_smpl          -   Number of model parameter samples to evaluate
%                           the objective function
%       prior_smpl      -   Set of N_smpl samples obtained from the prior
%                           distribution of the model parameters
%       NmaxSns         -   Maximum number of sensors
%       NmaxAct         -   Maximum number of actuators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % X - Definition of the objective function

    NmaxSns=length(GridPtX);     % No sensors
    % N_smpl                     % No samples
    NmaxAct=length(GridPtX);     % No actuators
    Nvar=size(P_sns_eval,1);    % No variables

    % Generate the inequality linear constraints
    A=zeros(2*(NmaxAct+NmaxSns+1),(NmaxAct+NmaxSns+1));
    for i=1:1*(NmaxAct+NmaxSns+1)
        for j=1:(NmaxAct+NmaxSns+1)
            if i==j
            A(i,j)=-1;
            end
        end
    end
    for i=(NmaxAct+NmaxSns+2):2*(NmaxAct+NmaxSns+1)
        for j=1:(NmaxAct+NmaxSns+1)
            if (i-(NmaxAct+NmaxSns+1))==j
            A(i,j)=1;
            end
        end
    end
    A(:,end)=0;
    A(2*(NmaxAct+NmaxSns+1),:)=0;
    b = [0*ones(1,(NmaxAct+NmaxSns+1)),1*ones(1,(NmaxAct+NmaxSns+1))];
    b(end)=0;

    % Equality contraints  for the number of sensors
    Aeq = [ones(1,NmaxAct+NmaxSns) -1];
    beq = 0;

    % Lower and upper bounds
    lb = [];
    ub = [];

    % Nonlinear constraints
    nonlcon = [];

    % Initial solution
    candSol=[ones(NmaxAct,1)./NmaxAct.*(Nini/2); ...
        ones(NmaxSns,1)./NmaxSns.*(1*Nini/2); Nini];

    % Options
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display','iter',...
    'StepTolerance',1e-7,'FunctionTolerance',1e-7);

    % Objective function
    fun=@(x)objFun(x,NmaxSns,NmaxAct,N_smpl,Nvar,x_cost,...
        y_cost,P_sns_eval);

    % Minimization of the objective function considering the gradient (Jacobian)
    [x, fval] = fmincon(fun,candSol,A,b,Aeq,beq,lb,ub,nonlcon,options);
    save('./res/Opt_Sol','x','fval')

end