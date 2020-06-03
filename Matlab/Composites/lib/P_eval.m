% Function: Evaluation of P matrices
function P_sns_eval=P_eval(GridPtX,GridPtY,prior_smpl,N_smpl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input paramters:
%       GridPtX,GridPtY -   X and Y coordinates of the grid points of
%                           possible sensor/actuator locations
%       prior_smpl      -   Set of N_smpl samples obtained from the prior
%                           distribution of the model parameters
%       N_smpl          -   Number of model parameter samples to evaluate
%                           the objective function
%
%   Output parameters:
%       P_sns_eval      -   Model-gradient multiplication matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NmaxSns=length(GridPtX);
    NmaxAct=length(GridPtX);
    P_sns_eval=zeros(size(prior_smpl,2),size(prior_smpl,2),NmaxSns,NmaxAct,N_smpl);
    for l=1:NmaxAct % Actuators
        for j=1:NmaxSns % Sensors
            for i=1:N_smpl
                P_sns_eval(:,:,j,l,i)=P_fun_mat(GridPtX(j), GridPtY(j), GridPtX(l), ...
                          GridPtY(l), prior_smpl(i,1), prior_smpl(i,2), prior_smpl(i,3),...
                          prior_smpl(i,4));
            end
        end
    end
   save(fullfile(cd, './dat/P_sns_eval.mat'), 'P_sns_eval')
end