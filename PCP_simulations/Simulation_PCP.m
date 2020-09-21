% Simulation_PCP
%
% simple script that calls generate_SNtrials to obtain sets of standard +
% outlier trials using various parameters for erp model and noise - it
% then computes various measures to evaluate the performance of the PCP
% method, in particlar it calls limo_pcout testing which data are seen as
% outliers and evaluate how much weigthing changes the data.
%
% Classification
% %%%%%%%%%%%%%%%
%
%                       | good trial generated |  bad trial generated |
%                       | ------------------------------------------- |
%  good trial detected  |    true positive    |    false positive     |
%                       | ------------------------------------------- |
%  bad trial detected   |    false negative   |     true negative     |
%                       | ------------------------------------------- |
%
%                                                 TP*TN - FP*FN
%  Matthew Correlation Coefficient: = ------------------------------------                                  
%                                     sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
%
% Robustness
% %%%%%%%%%%%
% influence: the Pearson correlation relative to the clean ERP 
% Kolmogorov-Smirnov test and distance: distance between the true average
%                                       and the one with outliers
%
% Cyril Pernet August 2019
% noise mixing developeed by Ignocio Muas
% ----------------------------------------

clear variables
nb_of_trials   = 200;
nb_of_outliers = [0.1*nb_of_trials 0.2*nb_of_trials 0.3*nb_of_trials ...
    0.4*nb_of_trials 0.5*nb_of_trials];
SNR = [1 1.5];

%% noise types

% after generating ERP on top of background, n percent are made outliers by
% adding to them a type of noise

noise{1} = 'white';
noise{2} = 'pink';
noise{3} = 'alpha';
noise{4} = 'gamma';

for n=1:4 % 4 types of noise
    
    MCC        = NaN(2,5,1000);
    ROC        = NaN(2,5,1000,2);
    influence1 = NaN(2,5,1000);
    influence2 = NaN(2,5,1000);
    influence3 = NaN(2,5,1000);
    diff1      = NaN(2,5,1000);
    ksstat1    = NaN(2,5,1000);
    diff2      = NaN(2,5,1000);
    ksstat2    = NaN(2,5,1000);
    diff3      = NaN(2,5,1000);
    ksstat3    = NaN(2,5,1000);
    
    for s = 1:2 % SNR
        for o = 1:5 % percentage of outliers
            fprintf('running %s noise: SNR %g %g%% outliers\n',noise{n},SNR(s),nb_of_outliers(o)/200*100);
            parfor MC = 1:1000
                
                [good_trials,outliers] = generate_SNtrials(nb_of_trials,nb_of_outliers(o), noise{n}, SNR(s), 1, 'off');
                good_trials            = good_trials';
                outliers               = outliers';
                all_trials             = [good_trials ; outliers];
                [~,out]                = limo_pcout(all_trials);
                good                   = [1:(nb_of_trials-nb_of_outliers(o))]';
                bad                    = [(nb_of_trials-nb_of_outliers(o))+1:200]';
                
                % classify
                TP           = length(intersect(good,find(out)));
                FP           = length(intersect(bad,find(out)));
                FN           = length(intersect(good,find(out==0)));
                TN           = length(intersect(bad,find(out==0)));
                if TP+FP+FN+TN ~= nb_of_trials
                    error('something is wrong in the classification')
                end
                ROC(s,o,MC,:)= [TP,FP]./nb_of_trials;
                MCC(s,o,MC)  = (TP*TN-FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
                
                % check robustness OLS/WLS/IRLS
                influence1(s,o,MC)                = corr(mean(all_trials)',mean(good_trials)');
                [diff1(s,o,MC),~,ksstat1(s,o,MC)] = kstest2(mean(all_trials),mean(good_trials));
                wavg                              = limo_WLS(ones(200,1),all_trials);
                influence2(s,o,MC)                = corr(wavg',mean(good_trials)');
                [diff2(s,o,MC),~,ksstat2(s,o,MC)] = kstest2(wavg,mean(good_trials));              
                wavg                              = limo_IRLS(ones(200,1),all_trials);
                influence3(s,o,MC)                = corr(wavg',mean(good_trials)');
                [diff3(s,o,MC),~,ksstat3(s,o,MC)] = kstest2(wavg,mean(good_trials));
            end
        end
    end
    
    if n==1
        save white_noise_sim
    elseif n == 2
        save pink_noise_sim
    elseif n == 3
        save alpha_noise_sim
    else
        save gamma_noise_sim
    end
end

    
%% amplitude of N1 noise

% after generating ERP on top of background, n percent are made outliers by
% making the N1 component smaller or bigger

OSR = [0.5 0.8 1.2 1.5];
disp('running N1 noise')
clear ROC MCC influence_curve difference KSSTAT

MCC        = NaN(2,5,4,1000);
ROC        = NaN(2,5,4,1000,2);
influence1 = NaN(2,5,4,1000);
influence2 = NaN(2,5,4,1000);
influence3 = NaN(2,5,4,1000);
diff1      = NaN(2,5,4,1000);
ksstat1    = NaN(2,5,4,1000);
diff2      = NaN(2,5,4,1000);
ksstat2    = NaN(2,5,4,1000);
diff3      = NaN(2,5,4,1000);
ksstat3    = NaN(2,5,4,1000);

for s = 1:2 % SNR
    for o = 1:5  % percentage of outliers 
        for sr = 1:4 % size of N1 outliers
            parfor MC = 1:1000
                [good_trials,outliers] = generate_SNtrials(nb_of_trials,nb_of_outliers(o), [], SNR(s), OSR(sr),'off');
                good_trials            = good_trials';
                outliers               = outliers';
                all_trials             = [good_trials ; outliers];
                [~,out]                = limo_pcout(all_trials);
                good                   = [1:(nb_of_trials-nb_of_outliers(o))]';
                bad                    = [(nb_of_trials-nb_of_outliers(o))+1:200]';
                
                % classify
                TP           = length(intersect(good,find(out)));
                FP           = length(intersect(bad,find(out)));
                FN           = length(intersect(good,find(out==0)));
                TN           = length(intersect(bad,find(out==0)));
                if TP+FP+FN+TN ~= nb_of_trials
                    error('something is wrong in the classification')
                end
                ROC(s,o,sr,MC,:)= [TP,FP]./nb_of_trials;
                MCC(s,o,sr,MC)  = (TP*TN-FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
                
                % check robustness OLS/WLS/IRLS
                influence1(s,o,sr,MC)                = corr(mean(all_trials)',mean(good_trials)');
                [diff1(s,o,sr,MC),~,ksstat1(s,o,MC)] = kstest2(mean(all_trials),mean(good_trials));
                wavg                                 = limo_WLS(ones(200,1),all_trials);
                influence2(s,o,sr,MC)                = corr(wavg',mean(good_trials)');
                [diff2(s,o,sr,MC),~,ksstat2(s,o,MC)] = kstest2(wavg,mean(good_trials));
                wavg                                 = limo_IRLS(ones(200,1),all_trials);
                influence3(s,o,sr,MC)                = corr(wavg',mean(good_trials)');
                [diff3(s,o,sr,MC),~,ksstat3(s,o,MC)] = kstest2(wavg,mean(good_trials));
            end
        end
    end
end
save amplitude_N1_noise_sim



