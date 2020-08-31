%% validation of repeated measures ANOVA

%% case 1 - 1 factor
% ------------------
% simulate 3 measures from a multivariate normal, with known means (2 3 6)
% and know variance structure such as a correlation exist between measures
% as expected from repeated measurements from the same subjects

N = 10; % number of subjects
mu = [2 3 6];
SIGMA = [1 0.5 0.2; 0.5 1 0.5; 0.2 0.5 1]; % data are correlated at .5 from 1 to 2 to 3 and .2 from 1 to 3
r = mvnrnd(mu,SIGMA,N); % get N observations per measures
xlswrite('1_factor',r)

data = nan(2,N,3); % data for limo 2 electrodes, N observations, 3 measures
data(1,:,:) = r;
data(2,:,:) = r;
gp = ones(N,1);
factors = [3];
result_Hotelling = limo_rep_anova(data,gp,factors) % standard Hotelling
result_robust_Hotelling = limo_robust_rep_anova(data,gp,factors,[],[],0) % robust Hotelling setting trimming at 0%

% same under H0, to check the type 1 error rate
P = NaN(1000,2);
for n=1:1000
    r = mvnrnd(mu,SIGMA,N);
    data = nan(2,N,3);
    data(1,:,:) = r;
    data(2,:,:) = r;
    gp = ones(N,1);
    factors = [3];
    result_Hotelling = limo_rep_anova(data,gp,factors);
    result_robust_Hotelling = limo_robust_rep_anova(data,gp,factors,[],[],0); % make sure to set percentage to 0 in limo_robust_rep_anova
    P(n,1) = result_Hotelling.p(1);
    P(n,2) = result_robust_Hotelling.p(1);
end

type1_error =  mean(P<0.05);


%% case 2 - N factor
% ------------------
% simulate 2*2 measures from a multivariate normal

mu = [2 3 6]; % same effect as above
mu = [mu mu*2]; % 1st factor is simply the data*2 and 2nd factor 
SIGMA = eye(6); 
S2 = [1 0.5 0.2; 0.5 1 0.5; 0.2 0.5 1]; % factor 2
SIGMA(1:3,1:3) = S2; SIGMA(4:6,4:6) = S2;

r = mvnrnd(mu,SIGMA,N); % get N observations per measures
xlswrite('2_factors',r)

data = nan(2,N,6); % data for limo_rep_anova
data(1,:,:) = r;
data(2,:,:) = r;
gp = ones(N,1);
factors = [2 3];
result2_Hotelling = limo_rep_anova(data,gp,factors) 
result2_robust_Hotelling = limo_robust_rep_anova(data,gp,factors,[],[],0) 

P = NaN(1000,6);
for n=1:1000
    r = mvnrnd(mu,SIGMA,N); % get N observations per measures
    data = nan(2,N,6); % data for limo_rep_anova
    data(1,:,:) = r;
    data(2,:,:) = r;
    gp = ones(N,1);
    factors = [2 3];
    result2_Hotelling = limo_rep_anova(data,gp,factors);
    result2_robust_Hotelling = limo_robust_rep_anova(data,gp,factors,[],[],0);
    P(n,1:3) = result2_Hotelling.p(:,1)';
    P(n,4:6) = result2_robust_Hotelling.p(:,1)';
end

type1_error =  mean(P<0.05);


% %% case 3 - gp*1 factor
% % ------------------
% % simulate 3 measures from a multivariate normal
% 
% mu = [2 3 6];
% SIGMA = [1 0.5 0.2; 0.5 1 0.5; 0.2 0.5 1]; % data are correlated at .5 from 1 to 2 to 3 and .2 from 1 to 3
% r1 = mvnrnd(mu,SIGMA,20); % get 20 observations per measures
% r2 = mvnrnd(mu,SIGMA,20);
% gp = [ones(20,1); ones(20,1)*2];
% r = [gp [r1 ; r2]]; xlswrite('gp_1_factors',r)
% 
% data = nan(2,40,3); % data for limo_rep_anova
% data(1,:,:) = r;
% data(2,:,:) = r;
% factors = [3];
% result = limo_rep_anova(data,gp,factors) 
% result0 = limo_robust_rep_anova(data,gp,factors) % make sure to set percentage to 0 in limo_robust_rep_anova
% 
% %% case 4 - gp*2 factors
% % ------------------
% % simulate 3 measures from a multivariate normal
% 
% mu = [2 3 6];
% SIGMA = [1 0.5 0.2; 0.5 1 0.5; 0.2 0.5 1]; % data are correlated at .5 from 1 to 2 to 3 and .2 from 1 to 3
% r1 = mvnrnd(mu,SIGMA,20); % get 20 observations per measures
% r2 = mvnrnd(mu,SIGMA,20);
% gp = [ones(20,1); ones(20,1)*2];
% r = [gp [r1 ; r2]]; xlswrite('gp_1_factors',r)
% 
% data = nan(2,40,3); % data for limo_rep_anova
% data(1,:,:) = r;
% data(2,:,:) = r;
% factors = [3];
% result = limo_rep_anova(data,gp,factors) 
% result0 = limo_robust_rep_anova(data,gp,factors) % make sure to set percentage to 0 in limo_robust_rep_anova
% 
