function plotOASP(X_d_v,Y_d_v,X_g_v,Y_g_v,prior_smpl,GridPtX,GridPtY,...
    minDist,NmaxSns,NmaxAct,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function that draws the optimal sensor and actuator configuration.
%
%   Input paramters:
%       X_d_v, Y_d_v    -   X and Y coordinates of the vertices of the
%                           area of possible damage occurrence
%       X_g_v, Y_g_v    -   X and Y coordinates of the vertices of the
%                           outer geometry of the plate
%       prior_smpl      -   Set of N_smpl samples obtained from the prior
%                           distribution of the model parameters
%       GridPtX,GridPtY -   X and Y coordinates of the grid points of
%                           possible sensor/actuator locations
%       minDist         -   Approximate distance between consecutive grid
%                           points (possible sensor/actuator locations)
%       NmaxSns         -   Maximum number of sensors
%       NmaxAct         -   Maximum number of actuators
%       x               -   Optimal solution containing z, w, and n
%
%   Output parameters:
%       This function only plots the results and save the figure in ./res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the resulting optimal layout
    figure;
    for i=1:length(X_d_v)-1
        oPol=plot([X_d_v(i),X_d_v(i+1)],[Y_d_v(i),Y_d_v(i+1)],...
             'color', [0.3,0.3,0.3]);
        hold on
        iPol=plot([X_g_v(i),X_g_v(i+1)],[Y_g_v(i),Y_g_v(i+1)],...
             'color', [0.,0.,0.]);
    end
    damSmpl=plot(prior_smpl(:,1),prior_smpl(:,2),'.','color', [0.7,0.7,0.7]);
    GridP=plot(GridPtX,GridPtY,'.');
    SnsVar=scatter(GridPtX,GridPtY,x(1:(NmaxSns))*100,'<','filled',...
        'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8,'MarkerFaceColor',[.1 .1 .8],...
        'MarkerEdgeColor',[.1 .1 .8]);
    ActVar=scatter(GridPtX,GridPtY,x((NmaxAct+1):end-1)*100,'>','filled',...
        'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6,'MarkerFaceColor',[.9 .7 .1],...
        'MarkerEdgeColor',[.9 .7 .1]);
    xlim([min(GridPtX)-minDist max(GridPtX)+minDist])
    ylim([min(GridPtY)-minDist max(GridPtY)+minDist])
    legend([SnsVar,ActVar],{'Sensors','Actuators'},...
        'interpreter','latex','fontsize',10,'location','northeast')
    ylabel('Y coordinate [m]','interpreter','latex','fontsize',10)
    xlabel('X coordinate [m]','interpreter','latex','fontsize',10)
    set(gca,'TickLabelInterpreter','latex','fontsize',10)
    set(gcf, 'Units', 'centimeters', 'OuterPosition', [5, 1, 15, 12]);
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'./res/Opt_Sol.pdf','-dpdf')
end