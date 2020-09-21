function [good_trials,outliers]=generate_SNtrials(nb_of_trials,nb_of_outliers, color, SNR, OSR,newfig)

% generate 0.5sec long trials (sampling 250Hz) with a P1 and N1 components
%
% FORMAT generate_trials(model,parameters,nb_of_trials,nb_of_outliers, color, OSR)
%
% INPUTS - nb of trials is the total number of trials to generate
%          (standards+outliers)
%        - nb of outliers is the total number of outliers to generate
%        - color is the type of noise to impose on outliers, if empty
%          outliers are defined based on amplitude of N1 otherwise they are
%          trials with addtional colored noise 'white' 'pink' 'alpha' 'theta'
%          'gamma' 'beta' 'linear trend'
%        - SNR is the average signal (event) to noise ratio
%        - OSR is the outlier to signal ratio (how much bigger outliers are
%          different from the other trials)
%        = newfig --> 'on' or 'off'
%
% OUTPUT good trials with amplitude based on input SNR
%        outlier trials with amplitude based on OSR (+color)
%
% The data generating process uses work from Yeung et al.
% standard trials have a power spectrum that matches human EEG
% REQUIRES https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator
%
% Cyril Pernet & Ignacio Suay Mas - August 2019
% -----------------------------------------------

%% check some inputs
if nb_of_outliers > nb_of_trials/2
    warndlg('there is more outliers than standard trials defined, from a distribution viewpoint they will then be the standards !')
end

if ~isempty(color) + ~strcmp(color,'white') + ~strcmp(color,'pink') ...
        + ~strcmp(color,'alpha') + ~strcmp(color,'theta') + ~strcmp(color,'gamma') ...
        + ~strcmp(color,'beta') + ~strcmp(color,'linear trend') == 8
    error('unknown noise color')
end

%% define defaults

% trials have 125 frames at 250Hz (0.5 sec) with a SNR of 0.2.
nb_events = nb_of_trials - nb_of_outliers;
frames    = 125;
srate     = 250;

%% simulate data

% background
background = reshape(noise(frames, nb_of_trials, srate),frames,nb_of_trials);

% SNR 
scaleP1 = (SNR*mean((max(background')))) / 1.5;
scaleN1 = (SNR*mean((min(background')))) * 1.5;

% events (P1, N1)
P1     = reshape(scaleP1.*peak(frames,nb_events,srate,15,35,2),frames,nb_events);
N1     = reshape(scaleN1.*peak(frames,nb_events,srate,15,45,2),frames,nb_events);
signal = (P1+N1);

% outliers
if isempty(color)
    scaleOutliers = OSR*mean(min(N1));
    P1            = reshape(scaleP1.*peak(frames,nb_of_outliers,srate,15,35,2),frames,nb_of_outliers);
    N1            = reshape((scaleOutliers).*peak(frames,nb_of_outliers,srate,15,45,2),frames,nb_of_outliers);
    outliers      = P1 + N1;
    
elseif strcmp(color,'white')
    P1                  = reshape(scaleP1.*peak(frames,nb_of_outliers,srate,15,35,2),frames,nb_of_outliers);
    N1                  = reshape(scaleN1.*peak(frames,nb_of_outliers,srate,15,45,2),frames,nb_of_outliers);
    hp                  = spectrum.periodogram('hamming');
    hpopts              = psdopts(hp,mean(N1)); set(hpopts,'Fs',srate);
    hpsd                = psd(hp,mean(N1,2),hpopts);
    power_freqdomain    = avgpower(hpsd);
    powerdb             = 10*log10(power_freqdomain/2);
    sigma               = abs(sqrt(powerdb / OSR));
    tmp                 = sigma*randn(frames,nb_of_outliers);
    noise_data          = zeros(frames,nb_of_outliers);
    index               = find(mean(P1+N1,2));
    index               = index(1):index(end)+length(index);
    noise_data(index,:) = tmp(index,:);
    outliers            = P1 + N1 + noise_data;
    
elseif strcmp(color,'pink')
    P1          = reshape(scaleP1.*peak(frames,nb_of_outliers,srate,15,35,2),frames,nb_of_outliers);
    N1          = reshape(scaleN1.*peak(frames,nb_of_outliers,srate,15,45,2),frames,nb_of_outliers);
    index       = find(mean(P1+N1,2)); 
    small_index = index(1):index(end);
    index       = index(1):index(end)+length(index);
    noise_data  = zeros(frames,nb_of_outliers);
    for j=1:length(index)
        [~, noise_data(index(j),:)] = phase_noise(nb_of_outliers,15,-SNR);
    end
    outliers    = P1 + N1 + noise_data;
    scale       = abs((mean(mean(outliers(small_index,:),2) ./ ...
        (mean(P1(small_index,:)+N1(small_index,:),2)))) / OSR);
    if scale > 1
        outliers = P1 + N1 + noise_data./scale;
    else
        outliers = P1 + N1 + noise_data.*scale;
    end
    
else % add some oscillations
    P1            = reshape(scaleP1.*peak(frames,nb_of_outliers,srate,15,35,2),frames,nb_of_outliers);
    N1            = reshape(scaleN1.*peak(frames,nb_of_outliers,srate,15,45,2),frames,nb_of_outliers);
    outliers      = P1 + N1;
    scaleOutliers = (OSR*range(mean(outliers,2)))*2;
    
    if strcmp(color,'theta');     f0=5;  
    elseif strcmp(color,'alpha'); f0=10;      
    elseif strcmp(color,'beta');  f0=20;       
    elseif strcmp(color,'gamma'); f0=40;         
    end
    
    index      = find(mean(P1+N1,2)); 
    index      = index(1):index(end)+length(index);
    time       = 1/srate:1/srate:0.5;
    noise_data = scaleOutliers*sin(2*pi*f0*time);  
    for i=1:nb_of_outliers
        delay = randi(round(frames/f0*2),1);   
        outliers(index,i) = (outliers(index,i)+noise_data(delay:delay+length(index)-1)')./2;
    end
end

good_trials = background(:,1:nb_events)     + signal;
outliers    = background(:,nb_events+1:end) + outliers;

% we can check it works by plotting the results
if strcmp(newfig,'on')
    figure; subplot(2,2,1); plot(1:125,background); title('background trials'); axis tight; grid on
    subplot(2,2,2); plot(1:125,good_trials); title('good trials'); axis tight; grid on
    subplot(2,2,3); plot(1:125,outliers); title('outliers'); axis tight; grid on
    subplot(2,2,4); plot(1:125,mean(good_trials,2),'LineWidth',2); hold on; axis tight; grid on
    boot_means = sort(squeeze(mean(reshape(good_trials(:,randi(nb_events,nb_events,600)),[125 nb_events 600]),2)),2);
    fillhandle = patch([[1:125] [125:-1:1]], [boot_means(:,15)' fliplr(boot_means(:,595)')], [0 0 1]);
    set(fillhandle,'EdgeColor',[0 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    plot(1:125,mean(outliers,2),'r','LineWidth',2); title('ERPs for good and outlier trials');
    boot_means = sort(squeeze(mean(reshape(outliers(:,randi(nb_of_outliers,nb_of_outliers,600)),[125 nb_of_outliers 600]),2)),2);
    fillhandle = patch([[1:125] [125:-1:1]], [boot_means(:,15)' fliplr(boot_means(:,595)')], [1 0 0]);
    set(fillhandle,'EdgeColor',[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
end

end