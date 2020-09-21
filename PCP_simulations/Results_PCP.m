% script to analyse the saved files from Simulation PCP
% generates plots for figures and summary stats csv files
% REQUIRES rst_data_plot https://github.com/CPernet/Robust_Statistical_Toolbox
%
% Cyril Pernet Septembre 2020
% ----------------------------

%% Classifcation Performences

% ROC space
% ----------------

figure('Name','ROC per % of colored outliers')
Percentage_noFP = NaN(4,2,5);
colour_map = limo_color_images(5);
for n = 1:4
    if n==1
        load white_noise_sim
    elseif n == 2
        load pink_noise_sim
    elseif n == 3
        load alpha_noise_sim
    else
        load gamma_noise_sim
    end
    
    for snr = 1:2
        if snr == 1
            subplot(2,4,n);
        else
            subplot(2,4,n+4);
        end
        plot(0:0.1:1,0:0.1:1,'k--')
        axis([0 1 0 1]); grid on; box on;
        hold on; xlabel('FP'); ylabel('TP');
        for o=1:5
            scatter(squeeze(ROC(snr,o,:,2)),squeeze(ROC(snr,o,:,1)),5,colour_map(o,:));
            Percentage_noFP(n,snr,o) = mean(squeeze(ROC(snr,o,:,2)) == 0)*100;
        end
        title([noise{n} ' noise: SNR = ' num2str(snr)])
    end
end

% check amplitude outliers
% -------------------------
load amplitude_N1_noise_sim
OSR = [0.5 0.8 1.2 1.5];
figure('Name','ROC per % of amplitude outliers')
for snr = 1:2
    for osr = 1:4
        if snr == 1
            subplot(2,4,osr);
        else
            subplot(2,4,osr+4);
        end
        plot(0:0.1:1,0:0.1:1,'k--')
        axis([0 1 0 1]); grid on; box on;
        hold on; xlabel('FP'); ylabel('TP');
        for o=1:5
            scatter(squeeze(ROC(snr,o,osr,:,2)),squeeze(ROC(snr,o,osr,:,1)),5,colour_map(o,:));
            Percentage_noFP(n,snr,o) = mean(squeeze(ROC(snr,o,:,2)) == 0)*100;
        end
        title(sprintf('SNR %g OSR %g',snr,OSR(osr)))
    end
end

% MCC analysis + plot TP/TN distributions
% ---------------------------------------
MedTP   = NaN(4,2,5);
TP_ci   = NaN(4,2,2,5);
MedFP   = NaN(4,2,5);
FP_ci   = NaN(4,2,2,5);
mcc_med = NaN(4,2,5);
mcc_ci  = NaN(4,2,2,5);

for n = 1:4
    figure('Name',[noise{n} ' noise']); 
   if n==1
        load white_noise_sim
    elseif n == 2
        load pink_noise_sim
    elseif n == 3
        load alpha_noise_sim
    else
        load gamma_noise_sim
   end
    
    plot_index = 1;
    for snr = 1:2
        subplot(2,3,plot_index);
        [MedTP(n,snr,:),TP_ci(n,snr,:,:)] = rst_data_plot(squeeze(ROC(snr,:,:,1))',...
            'estimator','median','bubble','off','point_size',10,'newfig','no');
        title([noise{n} ' noise: SNR = ' num2str(snr)]); ylabel('True Positive Rate'); drawnow
        subplot(2,3,plot_index+1);
        [MedFP(n,snr,:),FP_ci(n,snr,:,:)] = rst_data_plot(squeeze(ROC(snr,:,:,2))',...
            'estimator','median','bubble','off','point_size',10,'newfig','no');
        title([noise{n} ' noise: SNR = ' num2str(snr)]); ylabel('False Positive Rate'); drawnow
        subplot(2,3,plot_index+2);
        [mcc_med(n,snr,:),mcc_ci(n,snr,:,:)] = rst_data_plot(squeeze(MCC(snr,:,:))',...
            'estimator','median','bubble','off','point_size',10,'newfig','no');
        title([noise{n} ' noise: SNR = ' num2str(snr)]); ylabel('MCC'); drawnow
        plot_index = 4;
    end
end

% make a big table for these results
% ------------------------------
VN = {'True_Positive','False_Positive','MCC'};
RN = {'LSNR=1','HSNR=1','LSNR=2','HSNR=2'};
WN = table([squeeze(TP_ci(1,1,:,:)) ; squeeze(TP_ci(1,2,:,:))],...
    [squeeze(FP_ci(1,1,:,:)) ; squeeze(FP_ci(1,2,:,:))],...
    [squeeze(mcc_ci(1,1,:,:)) ; squeeze(mcc_ci(1,2,:,:))],'VariableNames',VN,'RowNames',RN);
disp('White noise 95% HDI'); disp(WN); writetable(WN,'PCP_WhiteNoise.csv')

PN = table([squeeze(TP_ci(2,1,:,:)) ; squeeze(TP_ci(2,2,:,:))],...
    [squeeze(FP_ci(2,1,:,:)) ; squeeze(FP_ci(2,2,:,:))],...
    [squeeze(mcc_ci(2,1,:,:)) ; squeeze(mcc_ci(2,2,:,:))],'VariableNames',VN,'RowNames',RN);
disp('Pink noise 95% HDI'); disp(PN); writetable(PN,'PCP_PinkNoise.csv')

AN = table([squeeze(TP_ci(3,1,:,:)) ; squeeze(TP_ci(3,2,:,:))],...
    [squeeze(FP_ci(3,1,:,:)) ; squeeze(FP_ci(3,2,:,:))],...
    [squeeze(mcc_ci(3,1,:,:)) ; squeeze(mcc_ci(3,2,:,:))],'VariableNames',VN,'RowNames',RN);
disp('Alpha noise 95% HDI'); disp(AN); writetable(AN,'PCP_AlphaNoise.csv')

GN = table([squeeze(TP_ci(4,1,:,:)) ; squeeze(TP_ci(4,2,:,:))],...
    [squeeze(FP_ci(4,1,:,:)) ; squeeze(FP_ci(4,2,:,:))],...
    [squeeze(mcc_ci(4,1,:,:)) ; squeeze(mcc_ci(4,2,:,:))],'VariableNames',VN,'RowNames',RN);
disp('Gamma noise 95% HDI'); disp(GN); writetable(GN,'PCP_GammaNoise.csv')

% check for amplitude
% -------------------
load amplitude_N1_noise_sim
for snr = 1:2
    figure; plot_index = 1;
    for n = 1:4
        subplot(4,3,plot_index);
        [MedTP(n,snr,:),TP_ci(n,snr,:,:)] = rst_data_plot(squeeze(ROC(snr,:,n,:,1))',...
            'estimator','median','bubble','off','point_size',10,'newfig','no');
        title(sprintf('SNR %g OSR %g',snr,OSR(n))); ylabel('True Positive Rate'); drawnow
        subplot(4,3,plot_index+1);
        [MedFP(n,snr,:),FP_ci(n,snr,:,:)] = rst_data_plot(squeeze(ROC(snr,:,n,:,2))',...
            'estimator','median','bubble','off','point_size',10,'newfig','no');
        title(sprintf('SNR %g OSR %g',snr,OSR(n))); ylabel('False Positive Rate'); drawnow
        subplot(4,3,plot_index+2);
        [mcc_med(n,snr,:),mcc_ci(n,snr,:,:)] = rst_data_plot(squeeze(MCC(snr,:,n,:))',...
            'estimator','median','bubble','off','point_size',10,'newfig','no');
        title(sprintf('SNR %g OSR %g',snr,OSR(n))); ylabel('MCC'); drawnow
        plot_index = plot_index + 3;
    end
end

WN = table([squeeze(TP_ci(1,1,:,:)) ; squeeze(TP_ci(1,2,:,:))],...
    [squeeze(FP_ci(1,1,:,:)) ; squeeze(FP_ci(1,2,:,:))],...
    [squeeze(mcc_ci(1,1,:,:)) ; squeeze(mcc_ci(1,2,:,:))],'VariableNames',VN,'RowNames',RN);
disp('0.5*N1 95% HDI'); disp(WN); writetable(WN,'PCP_05AmplitudeNoise.csv')

PN = table([squeeze(TP_ci(2,1,:,:)) ; squeeze(TP_ci(2,2,:,:))],...
    [squeeze(FP_ci(2,1,:,:)) ; squeeze(FP_ci(2,2,:,:))],...
    [squeeze(mcc_ci(2,1,:,:)) ; squeeze(mcc_ci(2,2,:,:))],'VariableNames',VN,'RowNames',RN);
disp('0.8*N1 95% HDI'); disp(PN); writetable(PN,'PCP_08AmplitudeNoise.csv')

AN = table([squeeze(TP_ci(3,1,:,:)) ; squeeze(TP_ci(3,2,:,:))],...
    [squeeze(FP_ci(3,1,:,:)) ; squeeze(FP_ci(3,2,:,:))],...
    [squeeze(mcc_ci(3,1,:,:)) ; squeeze(mcc_ci(3,2,:,:))],'VariableNames',VN,'RowNames',RN);
disp('1.2*N1 95% HDI'); disp(AN); writetable(AN,'PCP_12AmplitudeNoise.csv')

GN = table([squeeze(TP_ci(4,1,:,:)) ; squeeze(TP_ci(4,2,:,:))],...
    [squeeze(FP_ci(4,1,:,:)) ; squeeze(FP_ci(4,2,:,:))],...
    [squeeze(mcc_ci(4,1,:,:)) ; squeeze(mcc_ci(4,2,:,:))],'VariableNames',VN,'RowNames',RN);
disp('1.8*N1 95% HDI'); disp(GN); writetable(GN,'PCP_18AmplitudeNoise.csv')

%% Robustness

figure('Name','Pearson')
for n=1:4
    if n==1
        load white_noise_sim
    elseif n == 2
        load pink_noise_sim
    elseif n == 3
        load alpha_noise_sim
    else
        load gamma_noise_sim
    end
    
    for snr = 1:2
        
        [olsm, olsci]   = rst_pbCI(squeeze(influence1(snr,:,:))',599,0.5,'median');
        [wlsm, wlsci]   = rst_pbCI(squeeze(influence2(snr,:,:))',599,0.5,'median');
        [irlsm, irlsci] = rst_pbCI(squeeze(influence3(snr,:,:))',599,0.5,'median');
        
        if snr == 1
            subplot(2,4,n);
        else
            subplot(2,4,n+4);
        end
        title([noise{n} ' noise: snr ' num2str(snr)]); hold on; grid on; box on
        plot(olsm,'Color',colour_map(1,:)); fillhandle = patch([[1:5] [5:-1:1]], [olsci(1,:) fliplr(olsci(2,:))], colour_map(1,:));
        set(fillhandle,'EdgeColor',colour_map(1,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        plot(wlsm,'Color',colour_map(5,:)); fillhandle = patch([[1:5] [5:-1:1]], [wlsci(1,:) fliplr(wlsci(2,:))], colour_map(5,:));
        set(fillhandle,'EdgeColor',colour_map(5,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        plot(irlsm,'Color',colour_map(3,:)); fillhandle = patch([[1:5] [5:-1:1]], [irlsci(1,:) fliplr(irlsci(2,:))], colour_map(3,:));
        set(fillhandle,'EdgeColor',colour_map(3,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        axis([0.5 5.5 0.975 1])
    end
end

figure('Name','KS')
for n=1:4
    if n==1
        load white_noise_sim
    elseif n == 2
        load pink_noise_sim
    elseif n == 3
        load alpha_noise_sim
    else
        load gamma_noise_sim
    end
    
    for snr = 1:2
      
   [olsm, olsci]   = rst_pbCI(squeeze(ksstat1(snr,:,:))',599,0.5,'median');
   [wlsm, wlsci]   = rst_pbCI(squeeze(ksstat2(snr,:,:))',599,0.5,'median');
   [irlsm, irlsci] = rst_pbCI(squeeze(ksstat3(snr,:,:))',599,0.5,'median');
   
   if snr == 1
       subplot(2,4,n); 
   else
       subplot(2,4,n+4); 
   end
   title([noise{n} ' noise: snr ' num2str(snr)]); hold on; grid on; box on
        plot(olsm,'Color',colour_map(1,:)); fillhandle = patch([[1:5] [5:-1:1]], [olsci(1,:) fliplr(olsci(2,:))], colour_map(1,:));
        set(fillhandle,'EdgeColor',colour_map(1,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        plot(wlsm,'Color',colour_map(5,:)); fillhandle = patch([[1:5] [5:-1:1]], [wlsci(1,:) fliplr(wlsci(2,:))], colour_map(5,:));
        set(fillhandle,'EdgeColor',colour_map(5,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        plot(irlsm,'Color',colour_map(3,:)); fillhandle = patch([[1:5] [5:-1:1]], [irlsci(1,:) fliplr(irlsci(2,:))], colour_map(3,:));
        set(fillhandle,'EdgeColor',colour_map(3,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    end
end

load amplitude_N1_noise_sim
figure('Name','Pearson')
for n=1:4   
    for snr = 1:2
        
        [olsm, olsci]   = rst_pbCI(squeeze(influence1(snr,:,:))',599,0.5,'median');
        [wlsm, wlsci]   = rst_pbCI(squeeze(influence2(snr,:,:))',599,0.5,'median');
        [irlsm, irlsci] = rst_pbCI(squeeze(influence3(snr,:,:))',599,0.5,'median');
        
        if snr == 1
            subplot(2,4,n);
        else
            subplot(2,4,n+4);
        end
        title(sprintf('noise=%g osr=%g',snr,OSR(n)')); hold on; grid on; box on
        plot(olsm,'Color',colour_map(1,:)); fillhandle = patch([[1:5] [5:-1:1]], [olsci(1,:) fliplr(olsci(2,:))], colour_map(1,:));
        set(fillhandle,'EdgeColor',colour_map(1,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        plot(wlsm,'Color',colour_map(5,:)); fillhandle = patch([[1:5] [5:-1:1]], [wlsci(1,:) fliplr(wlsci(2,:))], colour_map(5,:));
        set(fillhandle,'EdgeColor',colour_map(5,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        plot(irlsm,'Color',colour_map(3,:)); fillhandle = patch([[1:5] [5:-1:1]], [irlsci(1,:) fliplr(irlsci(2,:))], colour_map(3,:));
        set(fillhandle,'EdgeColor',colour_map(3,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        axis([0.5 5.5 0.975 1])
    end
end

figure('Name','KS')
for n=1:4    
    for snr = 1:2
      
   [olsm, olsci]   = rst_pbCI(squeeze(ksstat1(snr,:,:))',599,0.5,'median');
   [wlsm, wlsci]   = rst_pbCI(squeeze(ksstat2(snr,:,:))',599,0.5,'median');
   [irlsm, irlsci] = rst_pbCI(squeeze(ksstat3(snr,:,:))',599,0.5,'median');
   
   if snr == 1
       subplot(2,4,n); 
   else
       subplot(2,4,n+4); 
   end
   title(sprintf('noise=%g osr=%g',snr,OSR(n)'));hold on; grid on; box on
        plot(olsm,'Color',colour_map(1,:)); fillhandle = patch([[1:5] [5:-1:1]], [olsci(1,:) fliplr(olsci(2,:))], colour_map(1,:));
        set(fillhandle,'EdgeColor',colour_map(1,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        plot(wlsm,'Color',colour_map(5,:)); fillhandle = patch([[1:5] [5:-1:1]], [wlsci(1,:) fliplr(wlsci(2,:))], colour_map(5,:));
        set(fillhandle,'EdgeColor',colour_map(5,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        plot(irlsm,'Color',colour_map(3,:)); fillhandle = patch([[1:5] [5:-1:1]], [irlsci(1,:) fliplr(irlsci(2,:))], colour_map(3,:));
        set(fillhandle,'EdgeColor',colour_map(3,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    end
end

