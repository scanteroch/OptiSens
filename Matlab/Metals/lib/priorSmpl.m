% Function: Generate prior samples for isotropic materials
function prior_smpl=priorSmpl(N_smpl,X_d_v,Y_d_v,V_mean,V_std)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input paramters:
%       N_smpl          -   Number of model parameter samples to evaluate
%                           the objective function
%       X_d_v, Y_d_v    -   X and Y coordinates of the vertices of the
%                           area of possible damage occurrence
%       V_mean, V_std   -   Mean and standard deviation of the prior PDF of
%                           the wave propagation velocity
%
%   Output parameters:
%       prior_smpl      -   Set of N_smpl samples obtained from the prior
%                           distribution of the model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prior_smpl=zeros(N_smpl,3);
    for i=1:N_smpl
        flagIsIn=0;
        while ~flagIsIn
            prior_smpl(i,1)=min(X_d_v)+(max(X_d_v)-min(X_d_v)).*rand(1);
            prior_smpl(i,2)=min(Y_d_v)+(max(Y_d_v)-min(Y_d_v)).*rand(1);
            flagIsIn = inpolygon(prior_smpl(i,1),prior_smpl(i,2),X_d_v,Y_d_v);
        end
    end
    prior_smpl(:,3) = normrnd(V_mean,V_std,[N_smpl,1]);
    save(fullfile(cd, './dat/prior_smpl.mat'),'prior_smpl','N_smpl')
end