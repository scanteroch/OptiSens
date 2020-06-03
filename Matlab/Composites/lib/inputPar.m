function [X_g_v,Y_g_v,X_d_v,Y_d_v,minDist,GridPtX,GridPtY,...
    N_smpl,prior_smpl,P_sns_eval,x_cost,y_cost,Nini]=inputPar(softnm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function that generates the necessary data for the optimization of
%   sensors and actuators. The user-dependent data will be introduced by
%   means of dialogue boxes.
%
%   Input paramters:
%       softnm - name of the software
%
%   Output parameters:
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Count the number of files in ./dat
    cont=0;
    my_file = './dat/prior_smpl.mat';
    if isfile(my_file)
        cont=cont+1;
    end
    my_file = './dat/GridPoints.mat';
    if isfile(my_file)
        cont=cont+1;
    end
    my_file = './dat/P_sns_eval.mat';
    if isfile(my_file)
        cont=cont+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If clause to select which case we are looking for:
    if cont==3
        prompt = {[newline,'\fontsize{9}','\it',sprintf(['HELP: The folder ./dat',...
            ' contains the data files from a previous study.\n\n']),...
            '\rmDo you want to reuse them? (Y/N)']};
            dlgtitle = [softnm,' - Initial point of fmincon'];
            definput = {'Y'};
            opts.Interpreter = 'tex';
            dlt = inputdlg(prompt,dlgtitle,[1, 70],definput,opts);
        while dlt{1}(1)~='Y' && dlt{1}(1)~='N'
            prompt = {[newline,'\fontsize{9}','\it',sprintf(['HELP: The folder ./dat',...
                ' contains the data files from a previous study.\n\n']),...
                '\rm\bf',sprintf('Note: Please, enter a valid answer.\n\n'),...
                '\rmDo you want to reuse them? (Y/N)']};
                dlgtitle = [softnm,' - Initial point of fmincon'];
                definput = {'Y'};
                opts.Interpreter = 'tex';
                dlt = inputdlg(prompt,dlgtitle,[1, 70],definput,opts);
        end
        if dlt{1}=='Y'  % Reuse files
            csIDX=1;
        else            % Delete files and start over
            csIDX=0;
        end
    elseif cont<3       % Delete files (if any) and start over
        csIDX=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch csIDX
        case 0          % Delete files
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Delete the files (if they exist)
            my_file = './dat/prior_smpl.mat';
            if isfile(my_file)
                delete(my_file)
            end
            my_file = './dat/GridPoints.mat';
            if isfile(my_file)
                delete(my_file)
            end
            my_file = './dat/P_sns_eval.mat';
            if isfile(my_file)
                delete(my_file)
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % X - Geometry of the plate
            prompt = {['\fontsize{9}','\it',sprintf('HELP: Both arrays need to have the same length.\n\n'),...
                '\rm 1. Enter an array of values for X coordinates in [m]:'],...
                '\fontsize{9} 2. Enter an array of values for Y coordinates in [m]:'};
            dlgtitle = [softnm,' - Outer geometry of the plate'];
            definput = {'[-1  ,  1  ,  1   , -0.9]','[-0.5, -0.5, -0.15,  0.5]'};
            opts.Interpreter = 'tex';
            X_Y_g_v = inputdlg(prompt,dlgtitle,[1, 70],definput,opts);
            X_g_v_aux=str2num(X_Y_g_v{1});
            Y_g_v_aux=str2num(X_Y_g_v{2});
            X_g_v=[X_g_v_aux X_g_v_aux(1)];
            Y_g_v=[Y_g_v_aux Y_g_v_aux(1)];
            clearvars X_Y_g_v X_g_v_aux Y_g_v_aux
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % X - Prior samples of the ToF model parameters (X_d,Y_d, and V)
            % X_g_v and X_d_v need to have the same number of vertices
            prompt = {[newline,'\fontsize{9}','\bf',sprintf(['First, enter information related to X_d and Y_d:\n\n']),...
                '\rm\it',sprintf(' HELP: Both arrays (1 and 2) need to have the same length between them, and the same number of vertices than the outer geometry.\n\n'),...
                '\rm 1. Enter an array of values for X coordinates in [m]:'],...
                '\fontsize{9} 2. Enter an array of values for Y coordinates in [m]:',...
                ['\fontsize{9}','\bf',sprintf('Second, enter the ellipse-based wave propagation velocity information:\n\n'),...
                '\rm',' 3. Enter the mean value of the horizontal axis in [m/s]'],...
                '\fontsize{9} 4. Enter the standard deviation of the horizontal axis:',...
                '\fontsize{9} 5. Enter the mean value of the vertical axis in [m/s]:',...
                '\fontsize{9} 6. Enter the standard deviation of the vertical axis:',...
                ['\fontsize{9}','\bf',...
                sprintf('Third, enter number of samples used for the objective function evaluation:\n\n'),...
                '\rm 7. Enter number of samples:']};
            dlgtitle = [softnm,' - Model parameters information'];
            definput = {'[-0.70, 0.70, 0.70 ,-0.65]','[-0.35,-0.35,-0.23 , 0.22]',...
                '6030','40','3549','40','10'};
            opts.Interpreter = 'tex';
            X_Y_d_v_V = inputdlg(prompt,dlgtitle,[1, 75],definput,opts);
            X_d_v_aux = str2num(X_Y_d_v_V{1});
            Y_d_v_aux = str2num(X_Y_d_v_V{2});
            a_mean = str2double(X_Y_d_v_V{3});
            a_std = str2double(X_Y_d_v_V{4});
            b_mean = str2double(X_Y_d_v_V{5});
            b_std = str2double(X_Y_d_v_V{6});
            N_smpl = str2double(X_Y_d_v_V{7}); % Number of samples for the MC integrals
            X_d_v=[X_d_v_aux X_d_v_aux(1)];
            Y_d_v=[Y_d_v_aux Y_d_v_aux(1)];
            clearvars X_Y_d_v_V X_d_v_aux Y_d_v_aux

            % Calculate or load the samples from the prior PDF
            my_file = './dat/prior_smpl.mat';
            if isfile(my_file)
                load(my_file);
            else
                prior_smpl=priorSmpl(N_smpl,X_d_v,Y_d_v,a_mean,a_std,...
                    b_mean, b_std);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % X - Grid of possible sensor locations
            prompt = {[newline,'\fontsize{9}','\it',sprintf(['HELP: These variables will create the grid in the',...
            ' space between the outer plate edges and the area of possible damage occurrence.\n\n']),...
                '\rm1. Enter the number of rings around the area of possible damage:'],...
                '\fontsize{9}2. Enter the minimum distance between consecutive points in the same ring in [m]:'};
            dlgtitle = [softnm,' - Grid of possible sensor locations'];
            definput = {'2','0.2'};
            opts.Interpreter = 'tex';
            Grid_Dist = inputdlg(prompt,dlgtitle,[1, 70],definput,opts);
            nSubGrid=str2double(Grid_Dist{1});
            minDist=str2double(Grid_Dist{2});
            clearvars Grid_Dist


            % Calculate or load the grid of possible sensor locations
            my_file = './dat/GridPoints.mat';
            if isfile(my_file)
                load(my_file)
            else
                [GridPtX, GridPtY]=gridPoints(nSubGrid,minDist,X_d_v,Y_d_v,X_g_v,Y_g_v);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Evaluation of the P matrices - Isotropic material

            % Calculate or load the grid of possible sensor locations
            my_file = './dat/P_sns_eval.mat';
            if isfile(my_file)
                load(my_file)
            else
                P_sns_eval=P_eval(GridPtX,GridPtY,prior_smpl,N_smpl);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Cost function definition by interpolating points (for Pchip) in [0,1]
            prompt = {[newline,'\fontsize{9}','\it',sprintf('HELP: Enter interpolating points used to build up the cost function with PCHIP function.\n\n'),...
                '\rm 1. Enter an array of values for X coordinates in number of sensors and actuators:'],...
                '\fontsize{9} 2. Enter an array of values for Y coordinates in dimensionless cost in the interval [0,1]:'};
            dlgtitle = [softnm,' - Definition of the cost function'];
            definput = {'[0  , 30 , 55, 60]','[0, 0.3, 0.95,1 ]'};
            opts.Interpreter = 'tex';
            costF = inputdlg(prompt,dlgtitle,[1, 70],definput,opts);
            x_cost=str2num(costF{1});
            y_cost=str2num(costF{2});
            clearvars costF
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Initial value
            prompt = {[newline,'\fontsize{9}','\it',sprintf(['HELP: It is recommended to use a higher number than',...
            ' the expected one.\n\n']),...
                '\rmEnter the initial guess of the number of sensors and actuators:']};
            dlgtitle = [softnm,' - Initial point of fmincon'];
            definput = {'30'};
            opts.Interpreter = 'tex';
            Nini = inputdlg(prompt,dlgtitle,[1, 70],definput,opts);
            Nini = str2double(Nini{1});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 1          % Reuse files
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load the files
            my_file = './dat/prior_smpl.mat';
            load(my_file)
            my_file = './dat/GridPoints.mat';
            load(my_file)
            my_file = './dat/P_sns_eval.mat';
            load(my_file)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Cost function definition by interpolating points (for Pchip) in [0,1]
            prompt = {[newline,'\fontsize{9}','\it',sprintf('HELP: Enter interpolating points used to build up the cost function with PCHIP function.\n\n'),...
                '\rm 1. Enter an array of values for X coordinates in number of sensors and actuators:'],...
                '\fontsize{9} 2. Enter an array of values for Y coordinates in dimensionless cost in the interval [0,1]:'};
            dlgtitle = [softnm,' - Definition of the cost function'];
            definput = {'[0  , 30 , 55, 60]','[0, 0.3, 0.95,1 ]'};
            opts.Interpreter = 'tex';
            costF = inputdlg(prompt,dlgtitle,[1, 70],definput,opts);
            x_cost=str2num(costF{1});
            y_cost=str2num(costF{2});
            clearvars costF
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Initial value
            prompt = {[newline,'\fontsize{9}','\it',sprintf(['HELP: It is recommended to use a higher number than',...
            ' the expected one.\n\n']),...
                '\rmEnter the initial guess of the number of sensors and actuators:']};
            dlgtitle = [softnm,' - Initial point of fmincon'];
            definput = {'30'};
            opts.Interpreter = 'tex';
            Nini = inputdlg(prompt,dlgtitle,[1, 70],definput,opts);
            Nini = str2double(Nini{1});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end