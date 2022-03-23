%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \sigma-k experimental surface maps for 
% STIRAP and saSTIRAP. 
% STIRAP: level 0,1,2 populations:
% p0_exp_STIRAP, p1_exp_STIRAP, p2_exp_STIRAP
% saSTIRAP: level 0,1,2 populations:
% p0_exp_saSTIRAP, p1_exp_saSTIRAP, p2_exp_saSTIRAP
%%%%%%%% Data from the experiments %%%%%%%%%%%%%
load('Exp_STIRAP.mat')
load('Exp_saSTIRAP.mat')
%%%%%%%%%%%%%%%%%%%%%%%

figure
 surf(sg*1e9, k, p2_exp_STIRAP); shading flat; view([-105 30]);
    xlabel('\sigma (ns)');
    ylabel('k');
    title('State |2\rangle population');
    xlim([sg(1)*1e9 sg(end)*1e9]);
    ylim([k(1) k(end)]);
    zlim([0,1])
    caxis([0 1])
    colorbar;
    
    hold on
    surf(sg*1e9, k, p2_exp_saSTIRAP); shading flat; view([-105 30]);
    xlabel('\sigma (ns)');
    ylabel('k');
    title('State |2\rangle population');
    xlim([sg(1)*1e9 sg(end)*1e9]);
    ylim([k(1) k(end)]);
    zlim([0,1])
    caxis([0 1])
    colorbar;

  


