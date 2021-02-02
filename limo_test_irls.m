function test = limo_test_irls(null)

%% IRLS validation
% ------------------
% this script imnplements a series of simple Monte Carlo simulation testing
% that the IRLS weights data as it should comparing with robustfit (this is
% not exactly the same implementation but both use a bisquare function and
% the same tuning contant, so this must be very similar.
% 
% FORMAT test = limo_test_irls(null)
% 
% INPUT  null is 'on' or 'off' (default) also testing the boostrap
% OUTPUT test is a matrix of RMSE difference to robustfit

NMC          = 1000;
sample_sizes = 10:10:100;
warning off % because robustfit will moan at rank deficiency
if nargin == 0
    null = 'off'; % default doesn't test boostrap undr H0 (ie bootstrapping NMC)
end
bootnull    = 1500; % if null is 'on' use 1000 boostraps

%% regression
MRMSE   = NaN(length(sample_sizes),2);
MCerror = NaN(length(sample_sizes),2);

parfor N = 1:length(sample_sizes)
    RMSE  = NaN(2,NMC);
    error = NaN(2,NMC);
    for MC = 1:NMC
        r = mvnrnd([0 0],eye(2),sample_sizes(N));
        Y = r(:,1);
        X = [r(:,2) ones(sample_sizes(N),1)];
        
        % if we use limo_glm with IRLS
        % ----------------------------
        [b,w] = limo_IRLS(X,Y);
        RMSE(1,MC) = sqrt(mean((Y-X*b).^2));
        if strcmpi(null,'on')
            Yhat        = X*b;
            Res         = Y-Yhat;
            HM          = (X.*w)*pinv(X.*w);
            df          = trace(HM'*HM)^2/trace((HM'*HM)*(HM'*HM))-1;
            dfe         = trace((eye(size(HM))-HM)'*(eye(size(HM))-HM));
            MSeffect    = norm((Yhat-mean(Yhat)))^2 / df;
            MSerror     = norm((Res-mean(Res)))^2 / dfe;
            F           = MSeffect / MSerror;
            p           = 1-fcdf(F,df,dfe);
            error(1,MC) = p <= 0.05;
            % null boostrap of this
            resampleb   = randi(sample_sizes(N),sample_sizes(N),bootnull);
            Yboot       = Y(resampleb);
            % WLS for X only
            B           = pinv(X.*w)*Yboot;
            Yhat        = X*B;
            Res         = Y-Yhat;
            MSeffect    = diag((Yhat-repmat(mean(Yhat),size(Y,1),1))'*(Yhat-repmat(mean(Yhat),size(Y,1),1)))./df;
            MSerror     = diag((Res-repmat(mean(Res),size(Y,1),1))'*(Res-repmat(mean(Res),size(Y,1),1)))./dfe;
            Fboot       = sort(MSeffect ./ MSerror);
            error(2,MC) = F >= Fboot(round(bootnull-5/100*bootnull));
        end
        
        % compare to robustfit
        % --------------------
        [~,stats] = robustfit(X(:,1) ,Y,'bisquare',4.685);
        RMSE(2,MC) = sqrt(mean(stats.resid.^2));
    end
    MRMSE(N,:)   = mean(RMSE,2);
    MCerror(N,:) = mean(error,2)';
end

figure('Name','IRLS test')
if strcmpi(null,'on')
    subplot(2,3,4);
    plot(sample_sizes,MCerror,'LineWidth',2);
    title(sprintf('Regression design\n %g Monte Carlo per samples sizes',NMC))
    grid on; box on; xlabel('Sample Size'); ylabel('error');
    legend('theoretical p','bootstrap F')
    subplot(2,3,1);
else
    subplot(1,3,1);
end
plot(sample_sizes,MRMSE,'LineWidth',2);
title(sprintf('Average RMSE - Regression design\n %g Monte Carlo per samples sizes',NMC))
grid on; box on; xlabel('Sample Size'); ylabel('RMSE'); legend('limo irls','robustfit')
test(:,1) = MRMSE(:,1)-MRMSE(:,2);


%% ANOVA/ANCOVA
MRMSE   = NaN(length(sample_sizes),4);
MCerror = NaN(length(sample_sizes),4);

parfor N = 1:length(sample_sizes)
    RMSE  = NaN(4,NMC);
    error = NaN(4,NMC);
    for MC = 1:1000
        r   = randn(sample_sizes(N)*3,4);
        Y   = r(1:sample_sizes(N),1:3);
        Y   = Y(:);
        cov = r(:,4);
        XX = [kron(eye(3),ones(sample_sizes(N),1)) cov ones(length(cov),1)];
        
        % if we use limo_glm with IRLS
        % ----------------------------
        X = [XX(:,[1 2 3]) XX(:,5)];
        [b,w] = limo_IRLS(X,Y);
        RMSE(1,MC) = sqrt(mean((Y-X*b).^2));
        if strcmpi(null,'on')
            Yhat        = X*b;
            Res         = Y-Yhat;
            HM          = (X.*w)*pinv(X.*w);
            df          = trace(HM'*HM)^2/trace((HM'*HM)*(HM'*HM))-1;
            dfe         = trace((eye(size(HM))-HM)'*(eye(size(HM))-HM));
            MSeffect    = norm((Yhat-mean(Yhat)))^2 / df;
            MSerror     = norm((Res-mean(Res)))^2 / dfe;
            F           = MSeffect / MSerror;
            p           = 1-fcdf(F,df,dfe);
            error(1,MC) = p <= 0.05;
            % null boostrap of this
            resampleb   = randi(sample_sizes(N)*3,sample_sizes(N)*3,bootnull);
            Yboot       = Y(resampleb);
            % WLS for X only
            B           = pinv(X.*w)*Yboot;
            Yhat        = X*B;
            Res         = Y-Yhat;
            MSeffect    = diag((Yhat-repmat(mean(Yhat),size(Y,1),1))'*(Yhat-repmat(mean(Yhat),size(Y,1),1)))./df;
            MSerror     = diag((Res-repmat(mean(Res),size(Y,1),1))'*(Res-repmat(mean(Res),size(Y,1),1)))./dfe;
            Fboot       = sort(MSeffect ./ MSerror);
            error(2,MC) = (1-mean(F>Fboot(round(bootnull-5/100*bootnull))))<=0.05;
        end
        
        X = XX; [b,w] = limo_IRLS(X,Y);
        RMSE(3,MC) = sqrt(mean((Y-X*b).^2));
        if strcmpi(null,'on')
            Yhat        = X*b;
            Res         = Y-Yhat;
            HM          = (X.*w)*pinv(X.*w);
            df          = trace(HM'*HM)^2/trace((HM'*HM)*(HM'*HM))-1;
            dfe         = trace((eye(size(HM))-HM)'*(eye(size(HM))-HM));
            MSeffect    = norm((Yhat-mean(Yhat)))^2 / df;
            MSerror     = norm((Res-mean(Res)))^2 / dfe;
            F           = MSeffect / MSerror;
            p           = 1-fcdf(F,df,dfe);
            error(3,MC) = p <= 0.05;
            % null boostrap of this
            resampleb   = randi(sample_sizes(N)*3,sample_sizes(N)*3,bootnull);
            Yboot       = Y(resampleb);
            % WLS for X only
            B           = pinv(X.*w)*Yboot;
            Yhat        = X*B;
            Res         = Y-Yhat;
            MSeffect    = diag((Yhat-repmat(mean(Yhat),size(Y,1),1))'*(Yhat-repmat(mean(Yhat),size(Y,1),1)))./df;
            MSerror     = diag((Res-repmat(mean(Res),size(Y,1),1))'*(Res-repmat(mean(Res),size(Y,1),1)))./dfe;
            Fboot       = sort(MSeffect ./ MSerror);
            error(4,MC) = 1-(mean(F > Fboot(round(bootnull-5/100*bootnull)))<=0.05);
        end
        
        % compare to robustfit
        % --------------------
        [~,stats] = robustfit(X(:,1:3) ,Y,'bisquare',4.685);
        %Cat = X(:,1:2); Cat(find(X(:,3)),:)=-1;
        %[B,stats] = robustfit(Cat ,Y,'bisquare',4.685);
        RMSE(2,MC) = sqrt(mean(stats.resid.^2));
        [~,stats] = robustfit(X(:,1:4) ,Y,'bisquare',4.685);
        %[B,stats] = robustfit([Cat X(:,4)] ,Y,'bisquare',4.685);
        RMSE(4,MC) = sqrt(mean(stats.resid.^2));
    end
    MRMSE(N,:)   = mean(RMSE,2);
    MCerror(N,:) = mean(error,2);
end
warning on

if strcmpi(null,'on')
    subplot(2,3,5);
    plot(sample_sizes,MCerror(:,[1 2]),'LineWidth',2);
    title(sprintf('ANOVA\n %g Monte Carlo per samples sizes',NMC))
    grid on; box on; xlabel('Sample Size'); ylabel('error');
    legend('theoretical p','bootstrap F')
    subplot(2,3,2);
else
    subplot(1,3,2);
end
plot(sample_sizes,MRMSE(:,[1 2]),'LineWidth',2);
title(sprintf('Average RMSE - ANOVA design\n %g Monte Carlo per samples sizes',NMC))
grid on; box on; xlabel('Sample Size'); ylabel('RMSE'); legend('limo irls','robustfit')

if strcmpi(null,'on')
    subplot(2,3,6);
    plot(sample_sizes,MCerror(:,[3 4]),'LineWidth',2);
    title(sprintf('ANCOVA\n %g Monte Carlo per samples sizes',NMC))
    grid on; box on; xlabel('Sample Size'); ylabel('error');
    legend('theoretical p','bootstrap F')
    subplot(2,3,3);
else
    subplot(1,3,3);
end

plot(sample_sizes,MRMSE(:,[3 4]),'LineWidth',2);
title(sprintf('Average RMSE - ANCOVA design\n %g Monte Carlo per samples sizes',NMC))
grid on; box on; xlabel('Sample Size'); ylabel('RMSE'); legend('limo irls','robustfit')
test(:,2) = MRMSE(:,1)-MRMSE(:,2);
test(:,3) = MRMSE(:,3)-MRMSE(:,4);

