function [avg_err,ci_err,max_err,maxci_err,maxc_err,maxcic_err] = limo_test_glmboot(chanlocs,varargin)

% routine the take a H0_folder generated by LIMO tools and return the error
% rate for the different files present - also make a figure to see if there
% is a consistency in space/frequency/time in the effect observed,
% typically due to non stationarity
%
% FORMAT [avg_err,ci_err] = limo_test_glmboot(H0,'step_size',100,'alphavalue',0.05,'figure','on')
%
% INPUTS H0 is a cell array of H0 folders (1 or many)
%        'alphavalue' is the nominal level to test
%        'step_size' the distance between resamples to look at convergence rate
%        'Nboot' the number of boostraps to use for error rate (e.g. use 1000 boot out of max computed)
%        'MinSamp' the lowest number of samples used for the null (from MinSamp to number of bootstrap-Nboot)
%        'figure' gives the option to check the space/frequency/time of
%                 effects - note if multiple H0 are given, data must all be of
%                 the same dimension
%
% OUTPUTS avg_err is the average error rate over the resamples
%         ci_err is the binomial 95% confidence interal
%
% Cyril Pernet January 2021

% defaults
alphavalue  = 0.05;
figurevalue = 'on';
MinSamp     = 600;
stepvalue   = 200;

if ischar(chanlocs)
    LIMO = load(chanlocs);
end

if ~iscell(varargin{1})
    error('cellarray expected as input')
else
    H0 = varargin{1};
    for key = 2:length(varargin)
        if ischar(varargin{key})
            value = key + 1;
            if contains(varargin{key},'alpha','IgnoreCase',true)
                alphavalue = varargin{value};
            elseif contains(varargin{key},'figure','IgnoreCase',true)
                figurevalue = varargin{value};
            elseif contains(varargin{key},'step','IgnoreCase',true)
                stepvalue = varargin{value};
            elseif contains(varargin{key},'Nboot','IgnoreCase',true)
                Nboot = varargin{value};
            elseif contains(varargin{key},'MinSamp','IgnoreCase',true)
                MinSamp = varargin{value};
            end
        end
    end
end

% get errors
for folder = length(H0):-1:1
    
    fileindex = 1;
    content = dir(fullfile(H0{folder},'H0_*.mat'));
    if isempty(content)
        fprintf('%s is empty\n',H0{folder})
    else
        fprintf('loading %s \n',H0{folder})
        for c=1:size(content,1)
            if ~contains(content(c).name,'Betas') && ~contains(content(c).name,'tfce') ...
                    && ~contains(content(c).name,'R2')
                data = load(fullfile(content(c).folder,content(c).name));
                data = data.(cell2mat(fieldnames(data)));
                if ~exist('Nboot','var')
                    Nboot = size(data,ndims(data));
                end
                % compute the number of error and average error
                if ndims(data) == 4
                    data(isnan(data(:,1,1,1)),:,:,:)   = [];
                    null_size                          = size(data,4)-Nboot;
                    if null_size < MinSamp
                        error('less than %g samples available to estimate the null - this is too low',MinSamp)
                    end
                    fprintf('Null data estimated from %g indepedent samples\n',null_size)
                    bootmax_position                   = zeros(size(data,1,2));
                    for b=null_size:-1:1               % use N-Nboot data for the null
                        tmp                            = squeeze(data(:,:,1,b));
                        [bootmax(b),p]                 = max(tmp(:));
                        [channel,time]                 = ind2sub([size(data,1) size(data,2)],p);
                        bootmax_position(channel,time) = bootmax_position(channel,time) + 1;
                    end
                    max_loc{folder,fileindex} = bootmax_position;
                    
                    disp('getting clusters under H0 boot ...');
                    bootcluster_position = zeros(size(data,1,2));
                    parfor boot = 1:null_size
                        % 1st find the cluster, thresholding H0 pvalues <= threshold p
                        [posclusterslabelmat,nposclusters] = limo_findcluster((squeeze(data(:,:,2,boot)) <= alphavalue),LIMO.channeighbstructmat,2); %#ok<PFBNS>
                        
                        % 2nd compute the mass for each cluster
                        bootM_b = squeeze(data(:,:,1,boot)); %#ok<PFBNS>
                        if nposclusters~=0
                            tmp = zeros(1,nposclusters);
                            for C = 1:nposclusters
                                tmp(C) = sum(bootM_b(posclusterslabelmat==C)); % sum stat value in a cluster label
                            end
                            [boot_maxclustersum(boot),clusterlabel] = max(tmp(:)); % save max value only
                            bootcluster_position = bootcluster_position+(posclusterslabelmat == clusterlabel);
                        else
                            boot_maxclustersum(boot) = 0;
                        end
                    end
                    cluster_loc{folder,fileindex} = bootcluster_position;
                    
                    % errors
                    Ntests                 = prod(size(data,[1 2]));
                    samp                   = null_size:-stepvalue:MinSamp; % size of the null
                    err{folder,fileindex}  = sum(squeeze(data(:,:,end,null_size+1:end)) < alphavalue,3); % number of errors of Nboot samples
                    for s = length(samp):-1:1
                        % cell-wise
                        tmp  = sum(squeeze(data(:,:,end,(end-samp(s)+1):end)) < alphavalue,3); % use samp instead of Nboot
                        [avg_err{folder,fileindex}(s),ci_err{folder,fileindex}(:,s)] = binofit(sum(tmp(:)),samp(s)*Ntests);
                        
                        fprintf('getting FWER based on %g bootstraps - subject %g file %g ...\n',samp(s),folder,fileindex)
                        maxs_err = zeros(1,Nboot); maxcs_err = maxs_err;
                        for b=size(data,4):-1:(null_size+1) % for each Nboot sample
                            % max
                            tmpmax              = bootmax;
                            tmpmax(tmpmax==Inf) = [];
                            tmpmax              = tmpmax(randperm(length(tmpmax)));
                            sortmaxM            = sort(tmpmax(1:samp(s)));
                            U                   = round((1-alphavalue).*samp(s));
                            mask                = single(squeeze(data(:,:,1,b)) >= sortmaxM(U));
                            if sum(mask(:)); maxs_err(b) = 1; end % indicates at least one error
                            clear tmpmax mask
                            
                            % cluster
                            tmpmax              = boot_maxclustersum;
                            tmpmax(tmpmax==Inf) = [];
                            tmpmax              = tmpmax(randperm(length(tmpmax)));
                            mask = limo_cluster_test(squeeze(data(:,:,1,b)),squeeze(data(:,:,2,b)), ...
                                tmpmax(1:samp(s)),LIMO.channeighbstructmat,2,alphavalue);
                            if sum(mask(:)); maxcs_err(b) = 1; end % indicates at least one error
                            clear tmpmax mask
                        end
                        [max_err{folder,fileindex}(s),maxci_err{folder,fileindex}(:,s)]   = binofit(sum(maxs_err),Nboot);
                        [maxc_err{folder,fileindex}(s),maxcic_err{folder,fileindex}(:,s)] = binofit(sum(maxcs_err),Nboot);
                    end
                    
                elseif ndims(data) == 5
                    data(isnan(data(:,1,1,1,1)),:,:,:,:) = [];
                    Ntests                               = prod(size(data,[1 2 3]));
                    samp                                 = size(data,5):-stepvalue:200;
                    bootmax_position                     = zeros(size(data,1,2,3));
                    for b=size(data,5):-1:1
                        tmp                                 = squeeze(data(:,:,:,1,b));
                        [bootmax(b),p]                      = max(tmp(:));
                        [channel,freq,time]                 = ind2sub([size(data,1) size(data,2) size(data,3)],p);
                        bootmax_position(channel,freq,time) = bootmax_position(channel,freq,time) + 1;
                    end
                    max_loc{folder,fileindex} = bootmax_position;
                    
                    disp('getting clusters under H0 boot ...');
                    bootcluster_position = zeros(size(data,1,2,3));
                    parfor boot = 1:size(data,5)
                        % 1st find the cluster, thresholding H0 pvalues <= threshold p
                        [posclusterslabelmat,nposclusters] = limo_findcluster((squeeze(data(:,:,:,2,boot)) <= alphavalue),LIMO.channeighbstructmat,2); %#ok<PFBNS>
                        
                        % 2nd compute the mass for each cluster
                        bootM_b = squeeze(data(:,:,:,1,boot)); %#ok<PFBNS>
                        if nposclusters~=0
                            tmp = zeros(1,nposclusters);
                            for C = 1:nposclusters
                                tmp(C) = sum(bootM_b(posclusterslabelmat==C)); % sum stat value in a cluster label
                            end
                            [boot_maxclustersum(boot),clusterlabel] = max(tmp(:)); % save max value only
                            bootcluster_position = bootcluster_position+(posclusterslabelmat == clusterlabel);
                        else
                            boot_maxclustersum(boot) = 0;
                        end
                    end
                    cluster_loc{folder,fileindex} = bootcluster_position;
                    
                    % errors
                    err{folder,fileindex}  = sum(squeeze(data(:,:,:,end,:)) < alphavalue,3); % number of errors
                    for s = length(samp):-1:1
                        % cell-wise
                        tmp  = sum(squeeze(data(:,:,:,end,1:samp(s))) < alphavalue,4);
                        [avg_err{folder,fileindex}(s),ci_err{folder,fileindex}(:,s)] = binofit(sum(tmp(:)),samp(s)*Ntests);
                        maxs_err = zeros(1,size(data,5)); maxcs_err = maxs_err;
                        fprintf('getting FWER based on %g bootstraps - subject %g file %g ...\n',samp(s),folder,fileindex)
                        for b=Nboot:-1:1
                            % max
                            tmpmax              = boot_maxclustersum;
                            tmpmax(b)           = [];
                            tmpmax(tmpmax==Inf) = [];
                            tmpmax              = tmpmax(randperm(length(tmpmax)));
                            if samp(s) > length(tmpmax)
                                tmpmax          = tmpmax(1:length(tmpmax));
                            else
                                tmpmax          = tmpmax(1:samp(s));
                            end
                            sortmaxM        = sort(tmpmax);
                            nboot           = length(sortmaxM);
                            U               = round((1-alphavalue).*nboot);
                            max_th          = sortmaxM(U);
                            mask            = squeeze(data(:,:,:,1,b)) >= max_th;
                            if sum(mask(:)); maxs_err(b) = 1; end % indicates at least one error
                            % cluster
                            tmpmax              = boot_maxclustersum;
                            tmpmax(b)           = [];
                            tmpmax(tmpmax==Inf) = [];
                            tmpmax              = tmpmax(randperm(length(tmpmax)));
                            if samp(s) > length(tmpmax)
                                tmpmax          = tmpmax(1:length(tmpmax));
                            else
                                tmpmax          = tmpmax(1:samp(s));
                            end
                            mask = limo_cluster_test(squeeze(data(:,:,:,1,b)),squeeze(data(:,:,:,2,b)), ...
                                tmpmax,LIMO.channeighbstructmat,2,alphavalue);
                            if sum(mask(:)); maxcs_err(b) = 1; end % indicates at least one error
                        end
                        [max_err{folder,fileindex}(s),maxci_err{folder,fileindex}(:,s)]   = binofit(sum(maxs_err),size(data,5));
                        [maxc_err{folder,fileindex}(s),maxcic_err{folder,fileindex}(:,s)] = binofit(sum(maxcs_err),size(data,5));
                    end
                else
                    error('null data files are expected to be 4 or 5 dimensionals')
                end
                
                if folder == 1
                    content(c).name(strfind(content(c).name,'_')) = ' ';
                    filename{fileindex} = content(c).name; % assuming each H0 is the same
                end
                fileindex = fileindex+1;
            end
        end
    end
end

% make figure
if strcmpi(figurevalue,'on')
    if size(err,2) == 1
        figure('Name','Error bias and rate');
        for type = 1:3
            if type == 1
                subplot(3,3,1);
                tmp = err{1,1};
            elseif type == 2
                subplot(3,3,2);
                tmp = max_loc{1,1};
            else
                subplot(3,3,3);
                tmp = cluster_loc{1,1};
            end
            
            if ndims(tmp) ==2 %#ok<ISMAT>
                imagesc(tmp)
            else
                imagesc(squeeze(mean(tmp,2)))
            end
            
            if type == 1
                title('Cell-wise Error density bias');
                ylabel('channels')
            elseif type == 2
                title('Max value Error density bias');
            else
                title('Cluster max Error density bias');
            end
        end
        
        cc = limo_color_images(size(err,1));
        for type = 1:3
            subplot(3,3,3+type); hold on
            for folder = size(err,1):-1:1 % subjects
                if ~isempty(err{folder,1})
                    if type == 1
                        low  = avg_err{folder,1}(end)  - ci_err{folder,1}(1,end);
                        high = ci_err{folder,1}(2,end) - avg_err{folder,1}(end);
                        errorbar(folder,avg_err{folder,1}(end),low,high,'LineWidth',2,'Color',cc(folder,:))
                        fprintf('subject''s error [%g %g]\n',ci_err{folder,1}(1,end),ci_err{folder,1}(2,end))
                    elseif type == 2
                        low  = max_err{folder,1}(end)   - maxci_err{folder,1}(1,end);
                        high = maxci_err{folder,1}(2,end) - max_err{folder,1}(end);
                        errorbar(folder,max_err{folder,1}(end),low,high,'LineWidth',2,'Color',cc(folder,:))
                        fprintf('subject''s error [%g %g]\n',maxci_err{folder,1}(1,end),maxci_err{folder,1}(2,end))
                    else
                        low  = maxc_err{folder,1}(end)     - maxcic_err{folder,1}(1,end);
                        high = maxcic_err{folder,1}(2,end) - maxc_err{folder,1}(end);
                        errorbar(folder,maxc_err{folder,1}(end),low,high,'LineWidth',2,'Color',cc(folder,:))
                        fprintf('subject''s error [%g %g]\n',maxcic_err{folder,1}(1,end),maxcic_err{folder,1}(2,end))
                    end
                end
            end
            plot(0:size(err,1),repmat(alphavalue,1,size(err,1)+1),'k')
            [~,ci]=binofit(round(Ntests*alphavalue),Ntests);
            plot(0:size(err,1),repmat(ci(1),1,size(err,1)+1),'k--')
            plot(0:size(err,1),repmat(ci(2),size(err,1)+1),'k--')
            grid on
            if type == 1
                title('Average cell-wise Type 1 Error');
                ylabel('subject''s type 1 error')
            elseif type == 2
                title('Type 1 FWER with max value correction');
            else
                title('Type 1 FWER with cluster mass correction');
            end
            
            subplot(3,3,6+type);
            LOW = zeros(1,length(avg_err{folder,1}));
            HIGH = LOW;
            for folder = size(err,1):-1:1
                if ~isempty(err{folder,1})
                    if type == 1
                        plot(samp,avg_err{folder,1},'--','LineWidth',1);hold on;
                        LOW  = LOW  + ci_err{folder,1}(1,:);
                        HIGH = HIGH + ci_err{folder,1}(2,:);
                        ylabel('type 1 error')
                    elseif type == 2
                        plot(samp,max_err{folder,1},'--','LineWidth',1);hold on;
                        LOW  = LOW  + maxci_err{folder,1}(1,:);
                        HIGH = HIGH + maxci_err{folder,1}(2,:);
                    else
                        plot(samp,maxc_err{folder,1},'--','LineWidth',1);hold on;
                        LOW  = LOW + maxcic_err{folder,1}(1,:);
                        HIGH = HIGH + maxcic_err{folder,1}(2,:);
                    end
                end
            end
            LOW  = LOW ./ (size(err,1)*size(err,2));
            HIGH = HIGH ./ (size(err,1)*size(err,2));
            plot(samp,LOW,'k','LineWidth',2);
            plot(samp,HIGH,'k','LineWidth',2);
            grid on; title('average convergence rate')
            fprintf('average error [%g %g]\n',LOW(end),HIGH(end))
        end
        
    else % do multiple figures
        figure('Name','cell-wise error density');
        t = tiledlayout('flow');
        for c=1:size(err,2)
            tmp = err{1,c};
            for folder = 2:size(err,1)
                if ~isempty(err{folder,c})
                    tmp = tmp + err{folder,c};
                end
            end
            nexttile(t);
            if ndims(tmp) ==2 %#ok<ISMAT>
                imagesc(tmp)
            else
                imagesc(squeeze(mean(tmp,2)))
            end
            title(sprintf('%s',filename{c}(1:end-4)))
            ylabel('type 1 error')
        end
        
        figure('Name','type 1 error')
        cc = limo_color_images(size(err,1));
        for c=1:size(err,2)
            subplot(size(err,2)+1,1,c); hold on
            for folder = size(err,1):-1:1
                if ~isempty(err{folder,c})
                    low = avg_err{folder,c}(end) - ci_err{folder,c}(1,end);
                    high = ci_err{folder,c}(2,end) - avg_err{folder,c}(end);
                    errorbar(folder,avg_err{folder,c}(end),low,high,'LineWidth',2,'Color',cc(folder,:))
                end
            end
            plot(0:size(err,1),repmat(alphavalue,1,size(err,1)+1),'k')
            [~,ci]=binofit(samp(end)*Ntests*0.05,samp(end)*Ntests);
            plot(0:size(err,1),repmat(ci(1),1,size(err,1)+1),'k--')
            plot(0:size(err,1),repmat(ci(2),size(err,1)+1),'k--')
            title(sprintf('%s',filename{c}(1:end-4))); grid on
        end
        
        subplot(size(err,2)+1,1,size(err,2)+1)
        AVG = zeros(1,length(avg_err{folder,c}));
        LOW = AVG ; HIGH = AVG;
        for c=1:size(err,2)
            for folder = 1:size(err,1)
                if ~isempty(err{folder,c})
                    AVG = AVG + avg_err{folder,c};
                    LOW = LOW + ci_err{folder,c}(1,:);
                    HIGH = HIGH + ci_err{folder,c}(2,:);
                end
            end
        end
        AVG  = AVG ./ (size(err,1)*size(err,2));
        LOW  = LOW ./ (size(err,1)*size(err,2));
        HIGH = HIGH ./ (size(err,1)*size(err,2));
        plot(samp,AVG,'LineWidth',2);
        hold on; plot(samp,LOW,'k--','LineWidth',1);
        plot(samp,HIGH,'k--','LineWidth',1);
        grid on; title('average convergence rate')
        ylabel('type 1 error')
    end
end