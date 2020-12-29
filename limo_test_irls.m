%% IRLS validation
% ------------------
% this script imnplement a series of simple Monte Carlo simulation testing
% that the IRLS weights data as it should comparing with robustfit (this is
% not exactly the same implementation but both use a bisquare function and
% the same tuning contant, so this must be very similar.

NMC = 1000;
sample_sizes = 10:10:100;
warning off % because robustfit will moan at rank deficiency

%% regression
MRMSE = NaN(length(sample_sizes),2);
parfor N = 1:length(sample_sizes)
    RMSE = NaN(2,NMC);
    for MC = 1:NMC
        r = mvnrnd([0 0],eye(2),sample_sizes(N));
        Y = r(:,1);
        X = [r(:,2) ones(sample_sizes(N),1)];
        
        % if we use limo_glm with IRLS
        % ----------------------------
        b = limo_IRLS(X,Y);
        RMSE(1,MC) = sqrt(mean((Y-X*b).^2));
        
        % compare to robustfit
        % --------------------
        [B,stats] = robustfit(X(:,1) ,Y,'bisquare',4.685);
        RMSE(2,MC) = sqrt(mean(stats.resid.^2));
    end
    MRMSE(N,:) = mean(RMSE,2);
end

figure('Name','RMSE')
subplot(1,3,1);
plot(sample_sizes,MRMSE,'LineWidth',2);
title(sprintf('Average RMSE - Regression design\n %g Monte Carlo per samples sizes',NMC))
grid on; box on; xlabel('Sample Size'); ylabel('RMSE'); legend('limo irls','robustfit')

%% ANOVA/ANCOVA
MRMSE = NaN(length(sample_sizes),4);
parfor N = 1:length(sample_sizes)
    RMSE = NaN(4,NMC);
    for MC = 1:1000
        r   = randn(sample_sizes(N)*3,4);
        Y   = r(1:sample_sizes(N),1:3);
        Y   = Y(:);
        cov = r(:,4);
        X = [kron(eye(3),ones(sample_sizes(N),1)) cov ones(length(cov),1)];
        
        % if we use limo_glm with IRLS
        % ----------------------------
        b = limo_IRLS([X(:,[1 2 3]) X(:,5)],Y);
        RMSE(1,MC) = sqrt(mean((Y-[X(:,[1 2 3]) X(:,5)]*b).^2));
        b = limo_IRLS(X,Y);
        RMSE(3,MC) = sqrt(mean((Y-X*b).^2));
        
        % compare to robustfit
        % --------------------
        [B,stats] = robustfit(X(:,1:3) ,Y,'bisquare',4.685);
        %Cat = X(:,1:2); Cat(find(X(:,3)),:)=-1;
        %[B,stats] = robustfit(Cat ,Y,'bisquare',4.685);
        RMSE(2,MC) = sqrt(mean(stats.resid.^2));
        [B,stats] = robustfit(X(:,1:4) ,Y,'bisquare',4.685);
        %[B,stats] = robustfit([Cat X(:,4)] ,Y,'bisquare',4.685);
        RMSE(4,MC) = sqrt(mean(stats.resid.^2));
    end
    MRMSE(N,:) = mean(RMSE,2);
end
warning on

subplot(1,3,2);
plot(sample_sizes,MRMSE(:,[1 2]),'LineWidth',2);
title(sprintf('Average RMSE - ANOVA design\n %g Monte Carlo per samples sizes',NMC))
grid on; box on; xlabel('Sample Size'); ylabel('RMSE'); legend('limo irls','robustfit')

subplot(1,3,3);
plot(sample_sizes,MRMSE(:,[3 4]),'LineWidth',2);
title(sprintf('Average RMSE - ANCOVA design\n %g Monte Carlo per samples sizes',NMC))
grid on; box on; xlabel('Sample Size'); ylabel('RMSE'); legend('limo irls','robustfit')
