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


    %%
    % saSTIRAP p2 map as a function of 
    % STIRAP area A (obtained from var1) and
    % saSTIRAP area A_02 (obtained from var2)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('Exp_saSTIRAP_area_map.mat')

    var1 = amplitude_02_range.^2;
    var2 = amplitude_01_range;
    title1 = 'Amplitude 02';
    title2 = 'Amplitude 01 (V)';
    
    figure
    surf(var2,var1,squeeze(state(:,:,1,3))); shading flat; view([90 -90]);
    xlabel(title2);
    ylabel(title1);
    xlim([var2(1), var2(end)]);
    ylim([var1(1), var1(end)]);
    zlim([0,1])
    caxis([0 1])
    colormap('hot')
    colorbar;
  


