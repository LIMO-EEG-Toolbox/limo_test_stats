% Simulation test PCP when one have less trials than data frames
%
% test the PCP algorithm for rank deficient cov. checking how the svd performs
% or alternatively using downsampling.
%
% The set of simulations performes the analysis with data generated at
% 1000/500/250Hz with 1501/1001/999/501/499/251/249/126 trials (10 to 50% white noise outliers)
% Importantly, to compare exactly thr weighting scheme, the same trials are
% reused, i.e. generate 1000Hz and 1501 trials then sample fron this.
%
% To quantify the effect of dimensionality and downsampling, one computes
% - Matthew Correlation Coefficient as a mesure of accuracy (how different results are)
% - hamming distance between outlier detection vectors (concordance among binary classification)
% - Pearson correlation between weights (concordance among weighting schemes)
%
% Cyril Pernet October 2021
% ----------------------------------------

clear variables
sampling      = [1000 500 250];
nb_of_trials  = [1501 1001 999 751 501 499 376 251 249 126];
for t = length(nb_of_trials):-1:1 
    nb_of_outliers(t,:) = floor([0.1*nb_of_trials(t) 0.2*nb_of_trials(t) 0.3*nb_of_trials(t) ...
        0.4*nb_of_trials(t) 0.5*nb_of_trials(t)]);
end


% after generating ERP on top of background, n percent are made outliers by
% adding to them white noise

MCC               = NaN(size(nb_of_outliers,2),length(nb_of_trials),3,1000);
ROV               = NaN(size(nb_of_outliers,2),length(nb_of_trials),3,2,1000);
classificationsim = NaN(size(nb_of_outliers,2),length(nb_of_trials),2,1000);
weightsim         = NaN(size(nb_of_outliers,2),length(nb_of_trials),2,1000);
O                 = NaN(size(nb_of_outliers,2),min(nb_of_trials),3,1000);
W                 = NaN(size(nb_of_outliers,2),min(nb_of_trials),3,1000);
 
for o = 1:size(nb_of_outliers,2)
    frames    = 1000;
    srate     = 1000;
    parfor MC = 1:1000
    fprintf('running: %g%% outliers MC %g\n',nb_of_outliers(1,o)/1501,MC);
                      
        % background
        background = reshape(noise(frames, max(nb_of_trials(:)), srate),frames,max(nb_of_trials(:)));
        % SNR = 1
        scaleP1 = mean((max(background'))) / 1.5;
        scaleN1 = mean((min(background'))) * 1.5;
        % events (P1, N1)
        nb_events = max(nb_of_trials(:))-nb_of_outliers(1,o);
        P1     = reshape(scaleP1.*peak(frames,nb_events,srate,60,140,8),frames,nb_events);
        N1     = reshape(scaleN1.*peak(frames,nb_events,srate,60,180,8),frames,nb_events);
        signal = (P1+N1);
        % add outliers (white noise)
        P1                  = reshape(scaleP1.*peak(frames,nb_of_outliers(1,o),srate,60,140,8),frames,nb_of_outliers(1,o));
        N1                  = reshape(scaleN1.*peak(frames,nb_of_outliers(1,o),srate,60,180,8),frames,nb_of_outliers(1,o));
        hp                  = spectrum.periodogram('hamming');
        hpopts              = psdopts(hp,mean(N1)); set(hpopts,'Fs',srate);
        hpsd                = psd(hp,mean(N1,2),hpopts);
        power_freqdomain    = avgpower(hpsd);
        powerdb             = 10*log10(power_freqdomain/2);
        sigma               = abs(sqrt(powerdb));
        tmp                 = sigma*randn(frames,nb_of_outliers(1,o));
        noise_data          = zeros(frames,nb_of_outliers(1,o));
        index               = find(mean(P1+N1,2));
        index               = index(1):index(end)+length(index);
        noise_data(index,:) = tmp(index,:);
        outliers            = P1 + N1 + noise_data;
        good_trials         = (background(:,1:nb_events)     + signal)';
        outliers            = (background(:,nb_events+1:end) + outliers)';
        
        % slice within parfoor 
        tmpA = NaN(length(nb_of_trials),2);
        tmpB = NaN(length(nb_of_trials),2);
        tmpC = NaN(length(nb_of_trials),3);
        tmpD = NaN(length(nb_of_trials),3,2);
        for s = 1:length(nb_of_trials)
            nb_events           = nb_of_trials(s)-nb_of_outliers(s,o);
            all_trials          = [good_trials(1:nb_events,:) ; outliers(1:nb_of_outliers(s,o),:)];
            [dist1,out1]        = limo_pcout(all_trials); % 1000Hz
            all_trials          = downsample(all_trials',2)';
            [dist2,out2]        = limo_pcout(all_trials); % 500Hz
            all_trials          = downsample(all_trials',2)';
            [dist3,out3]        = limo_pcout(all_trials); % 250Hz
            
            tmpA(s,:)         = [1-pdist([out1 out2]','hamming') 1-pdist([out1 out3]','hamming')];
            tmpB(s,:)         = [1-pdist([dist1 dist2]','correlation') 1-pdist([dist1 dist3]','correlation')];
               
            if s == length(nb_of_trials)
                 O(o,:,:,MC) = [out1 out2 out3];
                 W(o,:,:,MC) = [dist1 dist2 dist3];
            end
            
            % classify
            good         = 1:nb_events;
            bad          = nb_events+1:size(all_trials,1);
            out          = [out1 out2 out3];
            for ss = 1:3
                TP           = length(intersect(good,find(out(:,ss))));
                FP           = length(intersect(bad,find(out(:,ss))));
                FN           = length(intersect(good,find(out(:,ss)==0)));
                TN           = length(intersect(bad,find(out(:,ss)==0)));
                if TP+FP+FN+TN ~= nb_of_trials(s)
                    error('something is wrong in the classification')
                end
                tmpC(s,ss)   = (TP*TN-FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
                tmpD(s,ss,:) = [TP,FP]./nb_of_trials(s);
            end
        end
        % unsliced variables
        classificationsim(o,:,:,MC)  = tmpA;
        weightsim(o,:,:,MC)          = tmpB;
        MCC(o,:,:,MC)                = tmpC;
        ROC(o,:,:,:,MC)              = tmpD;
    end
    save pcp_limits_results
end

% does sampling affects performance
% --> for 40% of white noise outliers, lower sampling rate did better
figure('Name','Average Performance per sampling rate')
prop_outliers = round(nb_of_outliers(1,:)./max(nb_of_trials(:))*100);
TP = mean(squeeze(ROC(:,:,:,1,:)),4);
FP = mean(squeeze(ROC(:,:,:,2,:)),4);
subplot(1,4,1); plot(prop_outliers,squeeze(mean(TP,2)),'LineWidth',2); grid on; title('True Positives')
subplot(1,4,2); plot(prop_outliers,squeeze(mean(FP,2)),'LineWidth',2); grid on; title('False Positives')
subplot(1,4,3); hold on;
for d=1:3
    A = squeeze(FP(:,:,d));
    B = squeeze(TP(:,:,d));
    scatter(A(:),B(:)); 
end
grid on; title('ROC'); axis([0 1 0 1])
plot([0 1],[0 1],'k--','linewidth',2)
subplot(1,4,4); plot(prop_outliers,squeeze(mean((mean(MCC,4)),2)),'LineWidth',2); grid on; title('MCC')

% What is the effect of rank deficiency p>n
% --> up to 40% outliers, less trials means smaller MCC ; at 40% low
% sampling has a break down in performance and at 50% high sampling rate is better
figure
cc = limo_color_images(3);
thousand        = mean(squeeze(MCC(:,[1 2 3 5],sampling==1000,:)),3);
fivehundred     = mean(squeeze(MCC(:,[4 5 6 8],sampling==1000,:)),3);
twohundredfifty = mean(squeeze(MCC(:,[7 8 9 10],sampling==1000,:)),3);
for o=1:5
    subplot(1,5,o);
    plot([75 99 101 150], thousand(o,:),'lineWidth',2,'Color',cc(1,:))
    hold on; plot([75 99 101 150], fivehundred(o,:),'lineWidth',2,'Color',cc(2,:))
    plot([75 99 101 150], twohundredfifty(o,:),'lineWidth',2,'Color',cc(3,:))
    title([num2str(o*10) ' % outliers']); grid on
    if o == 1
        ylabel('Matthew corr. coef.')
    end
    xlabel('trials/time frame ratio')
end

figure('Name','Effect of Downsampling'); cindex = 1;
cc = limo_color_images(size(nb_of_outliers,2)*2);
subplot(2,2,1); hold on;
for o=1:size(nb_of_outliers,2)
    plot(fliplr(nb_of_trials),mean(squeeze(classificationsim(o,:,1,:)),2),'lineWidth',2,'Color',cc(cindex,:))
    cindex = cindex+1;
end
grid on; title('Classification hamming distance 1000Hz to 500Hz')
subplot(2,2,3); hold on;
for o=1:size(nb_of_outliers,2)
    plot(fliplr(nb_of_trials),mean(squeeze(classificationsim(o,:,2,:)),2),'lineWidth',2,'Color',cc(cindex,:))
    cindex = cindex+1;
end
grid on; title('Classification hamming distance 1000Hz to 250Hz')

subplot(2,2,2); hold on; cindex = 1;
for o=1:size(nb_of_outliers,2)
    plot(fliplr(nb_of_trials),mean(squeeze(weightsim(o,:,1,:)),2),'lineWidth',2,'Color',cc(cindex,:))
    cindex = cindex+1;
end
grid on; title('Weights correlations')
subplot(2,2,4); hold on;
for o=1:size(nb_of_outliers,2)
    plot(fliplr(nb_of_trials),mean(squeeze(weightsim(o,:,2,:)),2),'lineWidth',2,'Color',cc(cindex,:))
    cindex = cindex+1;
end
grid on; title('Weights correlations')

% we can also check for the last 126 trials when is always under sampled
MW = mean(W,4);  corr(squeeze(MO(1,:,:)))





