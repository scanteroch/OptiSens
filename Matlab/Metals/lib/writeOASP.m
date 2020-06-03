function writeOASP(x,NmaxAct,NmaxSns,GridPtX,GridPtY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function that writes a text file with the optimal sensor and actuator
%   configuration.
%
%   Input paramters:
%       x               -   Optimal solution containing z, w, and n
%       NmaxAct         -   Maximum number of actuators
%       NmaxSns         -   Maximum number of sensors
%       GridPtX,GridPtY -   X and Y coordinates of the grid points of
%                           possible sensor/actuator locations
%
%   Output parameters:
%       This function only writes the results in a text file in ./res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Table for the output
    Delim=   {'-------', '-------', '-------', '-------', '-------'};
    tabOutp_={'PZT #  '; 'X      '; 'Y      '; 'z      '; 'w      '};
    tabOutp =strjust(tabOutp_,'center');
    n_opt=round(x(end),2);
    z=x(1:NmaxSns); % Vector for the sensors
    w=x(NmaxAct+1:2*NmaxAct); % Vector for the actuators
    maxRow=round(max([n_opt,length(z(z>.2)),length(w(w>.2))])/2+1);
    [~, iDx]=sort(-z);
    fileID = fopen('./res/Output.txt','w');
    fprintf(fileID,"ENTROPY-BASED CONVEX OPTIMIZATION - ISOTROPIC PLATE\nOptimal sensor and actuator solution\n");
    fprintf(fileID,'\n');
    for i=1:length(tabOutp)
        fprintf(fileID,'%s ',tabOutp{i});
    end
    fprintf(fileID,'\n');
    for i=1:length(Delim)
        fprintf(fileID,'%s ',Delim{i});
    end
    fprintf(fileID,'\n');
    for i=1:maxRow
        fprintf(fileID,'%7.0f %7.3f %7.3f %7.3f %7.3f\n',i,GridPtX(iDx(i)),...
            GridPtY(iDx(i)),z(iDx(i)),w(iDx(i)));
    end
    fprintf(fileID,'\n');
    fprintf(fileID,'NOMENCLATURE\n');
    fprintf(fileID,'X   - X coordinate of the sensor/actuator [m]\n');
    fprintf(fileID,'Y   - Y coordinate of the sensor/actuator [m]\n');
    fprintf(fileID,'z   - Continuous decision variable for sensors\n');
    fprintf(fileID,'w   - Continuous decision variable for actuators\n\n');
    fprintf(fileID,'--------------------------------------------------\n\n');
    fprintf(fileID,'\nThe optimal number of sensors and actuators is: %.2f',n_opt);
    fclose(fileID);

end