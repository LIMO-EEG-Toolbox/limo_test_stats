function error_rate= limo_test_wglm(varargin)

% validation for OLS/WLS/IRLS
% simualate data under the null and run the analysis many times
% to compute the type 1 error rate
%
% FORMAT error_rate = limo_test_wglm(option)
%
% INPUTS correspond to the different cases of simuluation to perform
%        - Gaussian: simple normal data, no autocorrelation, no outliers
%        - EEG_bkg: eeg look alike data with P1 and N1
%
% to do it all: error_rate = limo_test_wglm('EEG_bkg','Gaussian')
%
% Cyril Pernet January 2021

if nargin == 0
    help limo_test_wglm
    return
else
    for c=1:nargin
        [ols,ci_ols,wls,ci_wls,irls,ci_irls] = compute_glm(varargin{c});
        error_rate.(varargin{c}).ols  = [ols,ci_ols];
        error_rate.(varargin{c}).wls  = [wls,ci_wls];
        error_rate.(varargin{c}).irls = [irls,ci_irls];
    end
end
end

function [ols,ci_ols,wls,ci_wls,irls,ci_irls] = compute_glm(option)

% how many Monte Carlo runs
nMC    = 10000;

% p-values place holders
pols   = NaN(100,nMC);
pwls   = NaN(100,nMC);
pirls  = NaN(100,nMC);

% erro rate place holders
erols  = NaN(100,nMC);
erwls  = NaN(100,nMC);
erirls = NaN(100,nMC);

%% Regression
clear X Y
fprintf('%s: running Regression model\n',option)
parfor MC = 1:nMC
    X = [randn(120,1) ones(120,1)];  % design matrix
    if strcmpi(option,'Gaussian')
        Y                      = randn(120,100); % generate 120 trials for 100 frames
    elseif strcmpi(option,'EEG_bkg')
        [good_trials,~]        = generate_SNtrials(120,0, [], 1, 1,'off');
        Y                      = good_trials(1:100,:)';
    else
        error('unknown option'); 
    end
    model        = limo_glm(Y, X, 0, 0, 1, 'OLS', 'Time', [], []);
    pols(:,MC)   = model.continuous.p;
    erols(:,MC)  = model.continuous.p<=0.05;
    model        = limo_glm(Y, X, 0, 0, 1, 'WLS', 'Time', [], []);
    pwls(:,MC)   = model.continuous.p;
    erwls(:,MC)  = model.continuous.p<=0.05;
    model        = limo_glm(Y, X, 0, 0, 1, 'IRLS', 'Time', [], []);
    pirls(:,MC)  = model.continuous.p;
    erirls(:,MC) = model.continuous.p<=0.05;
end

make_fig(option,'Regression',pols,pwls,pirls,erols,erwls,erirls)

% check confidence intervals
[ols.regression,ci_ols.regression]  = binofit(sum(erols(:)),numel(erols));
[wls.regression,ci_wls.regression]  = binofit(sum(erwls(:)),numel(erwls));
[irls.regression,ci_irls.regression] = binofit(sum(erirls(:)),numel(erirls));
fprintf('%s Regression error rates OLS [%g %g] WLS [%g %g] IRLS [%g %g]\n',...
    option, ci_ols.regression(1),ci_ols.regression(2),ci_wls.regression(1),...
    ci_wls.regression(2),ci_irls.regression(1),ci_irls.regression(2))


%% 1-way ANOVA
clear X Y
X = [kron(eye(3),ones(40,1)) ones(120,1)];  % design matrix
fprintf('%s: running ANOVA model\n',option)
parfor MC = 1:nMC
    if strcmpi(option,'Gaussian')
        Y                      = randn(120,100); % generate 120 trials for 100 frames
    elseif strcmpi(option,'EEG_bkg')
        [good_trials,~]        = generate_SNtrials(120,0, [], 1, 1,'off');
        Y                      = good_trials(1:100,:)';
    else
        error('unknown option'); 
    end
    model        = limo_glm(Y, X, 3, 0, 0, 'OLS', 'Time', [], []);
    pols(:,MC)   = model.conditions.p;
    erols(:,MC)  = model.conditions.p<=0.05;
    model        = limo_glm(Y, X, 3, 0, 0, 'WLS', 'Time', [], []);
    pwls(:,MC)   = model.conditions.p;
    erwls(:,MC)  = model.conditions.p<=0.05;
    model        = limo_glm(Y, X, 3, 0, 0, 'IRLS', 'Time', [], []);
    pirls(:,MC)  = model.conditions.p;
    erirls(:,MC) = model.conditions.p<=0.05;
end

make_fig(option,'ANOVA',pols,pwls,pirls,erols,erwls,erirls)

% check confidence intervals
[ols.anova,ci_ols.anova]  = binofit(sum(erols(:)),numel(erols));
[wls.anova,ci_wls.anova]  = binofit(sum(erwls(:)),numel(erwls));
[irls.anova,ci_irls.anova] = binofit(sum(erirls(:)),numel(erirls));
fprintf('%s ANOVA error rates OLS [%g %g] WLS [%g %g] IRLS [%g %g]\n',...
    option, ci_ols.anova(1),ci_ols.anova(2),ci_wls.anova(1),ci_wls.anova(2),...
    ci_irls.anova(1),ci_irls.anova(2))


%% 1-way ANCOVA

% covariate place holder
pols2   = NaN(100,nMC);
pwls2   = NaN(100,nMC);
pirls2  = NaN(100,nMC);
erols2  = NaN(100,nMC);
erwls2  = NaN(100,nMC);
erirls2 = NaN(100,nMC);

clear X Y
X = [kron(eye(3),ones(40,1)) randn(120,1) ones(120,1)];
fprintf('%s: running ANCOVA model',option)
parfor MC = 1:nMC
    if strcmpi(option,'Gaussian')
        Y                      = randn(120,100); % generate 120 trials for 100 frames
    elseif strcmpi(option,'EEG_bkg')
        [good_trials,~]        = generate_SNtrials(120,0, [], 1, 1,'off');
        Y                      = good_trials(1:100,:)';
    else
        error('unknown option'); 
    end
    model         = limo_glm(Y, X, 3, 0, 1, 'OLS', 'Time', [], []);
    pols(:,MC)    = model.conditions.p;
    erols(:,MC)   = model.conditions.p<0.05;
    pols2(:,MC)   = model.continuous.p;
    erols2(:,MC)  = model.continuous.p<0.05;
    model         = limo_glm(Y, X, 3, 0, 1, 'WLS', 'Time', [], []);
    pwls(:,MC)    = model.conditions.p;
    erwls(:,MC)   = model.conditions.p<=0.05;
    pwls2(:,MC)   = model.continuous.p;
    erwls2(:,MC)  = model.continuous.p<=0.05;
    model         = limo_glm(Y, X, 3, 0, 1, 'IRLS', 'Time', [], []);
    pirls(:,MC)   = model.conditions.p;
    erirls(:,MC)  = model.conditions.p<=0.05;
    pirls2(:,MC)  = model.continuous.p;
    erirls2(:,MC) = model.continuous.p<=0.05;
end

make_fig(option,'ANCOVA condition effect',pols,pwls,pirls,erols,erwls,erirls)

% check confidence intervals
[ols.ancova,ci_ols.ancova]  = binofit(sum(erols(:)),numel(erols));
[wls.ancova,ci_wls.ancova]  = binofit(sum(erwls(:)),numel(erwls));
[irls.ancova,ci_irls.ancova] = binofit(sum(erirls(:)),numel(erirls));
fprintf('%s ANCOVA condition error rates OLS [%g %g] WLS [%g %g] IRLS [%g %g]\n',...
    option,ci_ols.ancova(1),ci_ols.ancova(2),ci_wls.ancova(1),ci_wls.ancova(2),...
    ci_irls.ancova(1),ci_irls.ancova(2))

make_fig(option,'ANCOVA covariate effect',pols2,pwls2,pirls2,erols2,erwls2,erirls2)

[ols.cov,ci_ols.cov]  = binofit(sum(erols2(:)),numel(erols2));
[wls.cov,ci_wls.cov]  = binofit(sum(erwls2(:)),numel(erwls2));
[irls.cov,ci_irls.cov] = binofit(sum(erirls2(:)),numel(erirls2));
fprintf('%s ANCOVA covariate error rates OLS [%g %g] WLS [%g %g] IRLS [%g %g]\n',...
    option, ci_ols.cov(1),ci_ols.cov(2),ci_wls.cov(1),ci_wls.cov(2),ci_irls.cov(1),ci_irls.cov(2))

end

function make_fig(option,test,pols,pwls,pirls,erols,erwls,erirls)

figure('Name',[option ': Monte Carlo under H0']);
mp = [median(pols,2)';median(pwls,2)';median(pirls,2)']; % median p values over time
subplot(1,3,1)
scatter(mp(1,:),mp(2,:),50,[0 0 1])
xlabel('OLS'); ylabel('WLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(median(mp(:)),0.5,'ko','LineWidth',3)
plot(0.5,median(mp(:)),'ko','LineWidth',3)
title(sprintf('%s median p-values OLS vs WLS \n OLS error:%g, WLS error %g',test,median(erols(:)),median(erwls(:))))
subplot(1,3,2)
scatter(mp(1,:),mp(3,:),50,[0 0 1])
xlabel('OLS'); ylabel('IRLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(median(mp(:)),0.5,'ko','LineWidth',3)
plot(0.5,median(mp(:)),'ko','LineWidth',3)
title(sprintf('%s median p-values OLS vs IRLS \n OLS error:%g, IRLS error %g',test,median(erols(:)),median(erirls(:))))
subplot(1,3,3)
scatter(mp(2,:),mp(3,:),50,[0 0 1])
xlabel('WLS'); ylabel('IRLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(median(mp(:)),0.5,'ko','LineWidth',3)
plot(0.5,median(mp(:)),'ko','LineWidth',3)
title(sprintf('%s median p-values WLS vs IRLS \n WLS error:%g, IRLS error %g',test,median(erwls(:)),median(erirls(:))))
drawnow

end