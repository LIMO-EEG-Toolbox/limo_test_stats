%% IRLS validation
% ------------------
% this script imnplement a series of simple Monte Carlo simulation testing
% that the IRLS weights data as it should comparing with robustfit (this is
% not exactly the same implementation but both use a bisquare function and
% the same tuning contant, so on average this must be very similar

NMC = 1000;
sample_sizes = 10:10:100;

%% regression
parfor N = 1:length(sample_sizes)
    RMSE = NaN(2,NMC);
    for MC = 1:NMC
        r = mvnrnd([0 0],eye(2),sample_sizes(N));
        Y = r(:,1);
        X = [r(:,2) ones(sample_sizes(N),1)];
        
        % if we use limo_glm with IRLS
        % ----------------------------
        b = limo_IRLS(X,Y);
        RMSE(1,MC) = sqrt(mean((Y(1,:)'-X*b).^2));
        
        % compare to robustfit
        % --------------------
        B = robustfit(X(:,1) ,Y,'bisquare',4.685);
        RMSE(2,MC) = sqrt(mean((Y(1,:)'-X*flipud(B)).^2));
    end
    MRMSE(N,:) = mean(RMSE,2);
end

figure('Name','RMSE')
subplot(1,2,1);
plot(sample_sizes,MRMSE,'LineWidth',2);
title(sprintf('Regression design\nAverage RMSE over %g Monte Carlo per samples sizes',NMC))
grid on; box on; xlabel('Sample Size'); ylabel('RMSE'); legend('limo irls','robustfit')

%% ANCOVA

warning off % because robustfit will moan at rank deficiency
parfor N = 1:length(sample_sizes)
    RMSE = NaN(2,NMC);
    for MC = 1:1000
        r   = randn(sample_sizes(N)*3,4);
        Y   = r(1:sample_sizes(N),1:3);
        Y   = Y(:);
        cov = r(:,4);
        X = [kron(eye(3),ones(sample_sizes(N),1)) cov ones(length(cov),1)];
        
        % if we use limo_glm with IRLS
        % ----------------------------
        b = limo_IRLS(X,Y);
        RMSE(1,MC) = sqrt(mean((Y(1,:)'-X*b).^2));
        
        % compare to robustfit
        % --------------------
        B = robustfit(X(:,1:4) ,Y,'bisquare',4.685);
        RMSE(2,MC) = sqrt(mean((Y(1,:)'-X*flipud(B)).^2));
    end
    MRMSE(N,:) = mean(RMSE,2);
end
warning on

subplot(1,2,2);
plot(sample_sizes,MRMSE,'LineWidth',2);
title(sprintf('ANCOVA design\nAverage RMSE over %g Monte Carlo per samples sizes',NMC))
grid on; box on; xlabel('Sample Size'); ylabel('RMSE'); legend('limo irls','robustfit')
