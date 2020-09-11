% validates error rates for OLS/WLS/IRLS

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

%% 1-way ANOVA
X = [kron(eye(3),ones(40,1)) ones(120,1)];  % design matrix
parfor MC = 1:nMC
    Y            = randn(120,100); % generate 120 trials for 100 frames
    model        = limo_glm(Y, X, 3, 0, 0, 'OLS', 'Time', [], []);
    pols(:,MC)   = model.p;
    erols(:,MC)  = model.p<0.05;
    model        = limo_glm(Y, X, 3, 0, 0, 'WLS', 'Time', [], []);
    pwls(:,MC)   = model.conditions.p;
    erwls(:,MC)  = model.p<0.05;
    model        = limo_glm(Y, X, 3, 0, 0, 'IRLS', 'Time', [], []);
    pirls(:,MC)  = model.conditions.p;
    erirls(:,MC) = model.p<0.05;
end

figure('Name','Monte Carlo under H0');
mp = [mean(pols);mean(pwls);mean(pirls)]; % average frames
subplot(1,3,1)
scatter(mp(1,:),mp(2,:),50,[0 0 1])
xlabel('OLS'); ylabel('WLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(0.5,mean(mp(1,:)),'ko','LineWidth',3) 
plot(0.5,mean(mp(2,:)),'ko','LineWidth',3)
title(sprintf('ANOVA mean p-values OLS vs WLS \n OLS error:%g, WLS error %g',mean(mean(erols,2)),mean(mean(erwls,2))))
subplot(1,3,2)
scatter(mp(1,:),mp(3,:),50,[0 0 1])
xlabel('OLS'); ylabel('IRLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(0.5,mean(mp(1,:)),'ko','LineWidth',3)
plot(0.5,mean(mp(3,:)),'ko','LineWidth',3)
title(sprintf('ANOVA mean p-values OLS vs IRLS \n OLS error:%g, IRLS error %g',mean(mean(erols,2)),mean(mean(erirls,2))))
subplot(1,3,3)
scatter(mp(2,:),mp(3,:),50,[0 0 1])
xlabel('WLS'); ylabel('IRLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(0.5,mean(mp(1,:)),'ko','LineWidth',3)
plot(0.5,mean(mp(3,:)),'ko','LineWidth',3)
title(sprintf('ANOVA mean p-values WLS vs IRLS \n WLS error:%g, IRLS error %g',mean(mean(erwls,2)),mean(mean(erirls,2))))

%% 1-way ANCOVA
pols2   = NaN(100,nMC); 
pwls2   = NaN(100,nMC);
pirls2  = NaN(100,nMC);
erols2  = NaN(100,nMC);
erwls2  = NaN(100,nMC);
erirls2 = NaN(100,nMC);

X = [kron(eye(3),ones(40,1)) randn(120,1) ones(120,1)];
parfor MC = 1:nMC
    Y             = randn(120,100); 
    model         = limo_glm(Y, X, 3, 0, 1, 'OLS', 'Time', [], []);
    pols(:,MC)    = model.conditions.p;
    erols(:,MC)   = model.conditions.p<0.05;
    pols2(:,MC)   = model.continuous.p;
    erols2(:,MC)  = model.continuous.p<0.05;
    model         = limo_glm(Y, X, 3, 0, 1, 'WLS', 'Time', [], []);
    pwls(:,MC)    = model.conditions.p;
    erwls(:,MC)   = model.conditions.p<0.05;
    pwls2(:,MC)   = model.continuous.p;
    erwls2(:,MC)  = model.continuous.p<0.05;
    model         = limo_glm(Y, X, 3, 0, 1, 'IRLS', 'Time', [], []);
    pirls(:,MC)   = model.conditions.p;
    erirls(:,MC)  = model.conditions.p<0.05;
    pirls2(:,MC)  = model.continuous.p;
    erirls2(:,MC) = model.continuous.p<0.05;
end

figure('Name','Monte Carlo under H0');
mp = [mean(pols);mean(pwls);mean(pirls)]; % average frames
subplot(2,3,1)
scatter(mp(1,:),mp(2,:),50,[0 0 1])
xlabel('OLS'); ylabel('WLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(0.5,mean(mp(1,:)),'ko','LineWidth',3)
plot(0.5,mean(mp(2,:)),'ko','LineWidth',3)
title(sprintf('ANCOVA mean p-values OLS vs WLS \n OLS error:%g, WLS error %g',mean(mean(erols,2)),mean(mean(erwls,2))))
subplot(2,3,2)
scatter(mp(1,:),mp(3,:),50,[0 0 1])
xlabel('OLS'); ylabel('IRLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(0.5,mean(mp(1,:)),'ko','LineWidth',3)
plot(0.5,mean(mp(3,:)),'ko','LineWidth',3)
title(sprintf('ANCOVA mean p-values OLS vs IRLS \n OLS error:%g, IRLS error %g',mean(mean(erols,2)),mean(mean(erirls,2))))
subplot(2,3,3)
scatter(mp(2,:),mp(3,:),50,[0 0 1])
xlabel('WLS'); ylabel('IRLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(0.5,mean(mp(1,:)),'ko','LineWidth',3)
plot(0.5,mean(mp(3,:)),'ko','LineWidth',3)
title(sprintf('ANCOVA mean p-values WLS vs IRLS \n WLS error:%g, IRLS error %g',mean(mean(erwls,2)),mean(mean(erirls,2))))

mp = [mean(pols2);mean(pwls2);mean(pirls2)]; % average frames
subplot(2,3,4)
scatter(mp(1,:),mp(2,:),50,[0 0 1])
xlabel('OLS'); ylabel('WLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(0.5,mean(mp(1,:)),'ko','LineWidth',3)
plot(0.5,mean(mp(2,:)),'ko','LineWidth',3)
title(sprintf('Reg mean p-values OLS vs WLS \n OLS error:%g, WLS error %g',mean(mean(erols2,2)),mean(mean(erwls2,2))))
subplot(2,3,5)
scatter(mp(1,:),mp(3,:),50,[0 0 1])
xlabel('OLS'); ylabel('IRLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(0.5,mean(mp(1,:)),'ko','LineWidth',3)
plot(0.5,mean(mp(3,:)),'ko','LineWidth',3)
title(sprintf('Reg mean p-values OLS vs IRLS \n OLS error:%g, IRLS error %g',mean(mean(erols2,2)),mean(mean(erirls2,2))))
subplot(2,3,6)
scatter(mp(2,:),mp(3,:),50,[0 0 1])
xlabel('WLS'); ylabel('IRLS')
axis([0.4 0.6 0.4 0.6])
hold on
plot([0.4 0.6],[0.4 0.6],'r','LineWidth',2);
grid on; box on
plot(0.5,mean(mp(1,:)),'ko','LineWidth',3)
plot(0.5,mean(mp(3,:)),'ko','LineWidth',3)
title(sprintf('Reg mean p-values WLS vs IRLS \n WLS error:%g, IRLS error %g',mean(mean(erwls2,2)),mean(mean(erirls2,2))))
