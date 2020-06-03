% Function: Generate prior samples for isotropic materials
function [GridPtX, GridPtY]=gridPoints(nSubGrid,minDist,X_d_v,...
    Y_d_v,X_g_v,Y_g_v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input paramters:
%       nSubGrid        -   Number of rings around the area of possible 
%                           damage occurrence for the grid of points
%       minDist         -   Approximate distance between consecutive grid
%                           points (possible sensor/actuator locations)
%       X_d_v, Y_d_v    -   X and Y coordinates of the vertices of the
%                           area of possible damage occurrence
%       X_g_v, Y_g_v    -   X and Y coordinates of the vertices of the
%                           outer geometry of the plate
%
%   Output parameters:
%       GridPtX,GridPtY -   X and Y coordinates of the grid points of
%                           possible sensor/actuator locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the corner points
    auxGrdPtX=zeros(nSubGrid,length(X_d_v)+1);
    auxGrdPtY=zeros(nSubGrid,length(X_d_v)+1);
    for i=1:nSubGrid
        for j=1:length(X_d_v)
            auxGrdPtX(i,j)=X_g_v(j)-(i)*(X_g_v(j)-X_d_v(j))/(nSubGrid+1);
            auxGrdPtY(i,j)=Y_g_v(j)-(i)*(Y_g_v(j)-Y_d_v(j))/(nSubGrid+1);
        end
        auxGrdPtX(i,j+1)=auxGrdPtX(i,1);
        auxGrdPtY(i,j+1)=auxGrdPtY(i,1);
    end
    
    % Generate the grid points with approx minDist between them:
    GridPtX=[];
    GridPtY=[];
    for j=1:nSubGrid
        for i=1:length(X_d_v)-1
            
            auxDist=sqrt((auxGrdPtX(j,i+1)-auxGrdPtX(j,i)).^2+...
                            (auxGrdPtY(j,i+1)-auxGrdPtY(j,i)).^2);
            auxNpt=round(auxDist/minDist);
            
            if auxGrdPtX(j,i)==auxGrdPtX(j,i+1)
                GrdPtXtemp=linspace(auxGrdPtX(j,i),auxGrdPtX(j,i+1),auxNpt);
            else
                GrdPtXtemp=auxGrdPtX(j,i):(auxGrdPtX(j,i+1)-...
                    auxGrdPtX(j,i))/auxNpt:auxGrdPtX(j,i+1)-...
                    (auxGrdPtX(j,i+1)-auxGrdPtX(j,i))/auxNpt;
                if isempty(GrdPtXtemp)
                    GrdPtXtemp=auxGrdPtX(j,i);
                end
            end
            
            if auxGrdPtY(j,i)==auxGrdPtY(j,i+1)
                GrdPtYtemp=linspace(auxGrdPtY(j,i),auxGrdPtY(j,i+1),auxNpt);
            else
                GrdPtYtemp=auxGrdPtY(j,i):(auxGrdPtY(j,i+1)-...
                    auxGrdPtY(j,i))/auxNpt:auxGrdPtY(j,i+1)-...
                    (auxGrdPtY(j,i+1)-auxGrdPtY(j,i))/auxNpt;
                if isempty(GrdPtYtemp)
                    GrdPtYtemp=auxGrdPtY(j,i);
                end
            end
            GridPtX=[GridPtX GrdPtXtemp];
            GridPtY=[GridPtY GrdPtYtemp];
        end
    end
    save(fullfile(cd, './dat/GridPoints.mat'), 'GridPtX', 'GridPtY',...
        'X_g_v','Y_g_v','X_d_v','Y_d_v','minDist')
end