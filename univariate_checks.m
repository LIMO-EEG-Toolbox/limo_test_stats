%% IRLS validation
% ------------------
% this script imnplement a series of simple Monte Carlo simulation testing
% that the bootstrap procedures implemented in LIMO EEG gives the right
% type one error rate

%% Simple regression
% we compare OLS to IRLS with/without bootstrap
% case 1: data are normally distributed 

clear
% N=20 corr = 0
% -----------------------
warning off
for MC = 1:10
    r = mvnrnd([0 0],eye(2),20);
    Y(1,:) = r(:,1);
    Y(2,:) = r(:,1);
    X = [r(:,2) ones(20,1)];
    nb_conditions = 0;
    nb_continuous = 1;
    nb_interactions = 0;
    
    % if we use OLS
    % --------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'OLS');
    F_obs =  model.continuous.F(1);
    p_value_OLS(MC) = model.continuous.p(1);
    
    % if we use OLS and bootstrap
    % --------------------------
    boot_table = randi(20,20,599);
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,'OLS',boot_table);
    for b=1:599; bootF(b) =  model.continuous.F{b}(1); end
    p_value_OLS_boot(MC) = 1 - (sum(F_obs>bootF) / 599);
    
    % if we use IRLS
    % ---------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'IRLS');
    F_obs =  model.continuous.F(1);
    p_value_IRLS(MC) = model.continuous.p(1);
    
    % if we use IRLS and bootstrap
    % ----------------------------
    weights = model.W';
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,'IRLS',boot_table);
    for b=1:599; bootF(b) =  model.continuous.F{b}(1); end
    p_value_IRLS_boot(MC) = 1 - (sum(F_obs>bootF) / 599);
end
save single_regression;    

%% multiple regression
% we compare OLS to IRLS with/without bootstrap
% case 1: data are normally distributed 
clear

% N=20 corr = 0
% -----------------------
for MC = 1:1000
    r = mvnrnd([0 0 0],eye(3),20);
    Y(1,:) = r(:,1);
    Y(2,:) = r(:,1);
    X = [r(:,[2 3]) ones(20,1)];
    nb_conditions = 0;
    nb_continuous = 2;
    nb_interactions = 0;
    
    % if we use OLS
    % --------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'OLS');
    F1_obs =   model.continuous.F(1,1);
    F2_obs =   model.continuous.F(1,2);
    p_value1_OLS(MC) = model.continuous.p(1,1);
    p_value2_OLS(MC) = model.continuous.p(1,2);

    % if we use OLS and bootstrap
    % --------------------------
    weights = model.W';
    boot_table = randi(20,20,599);
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,weights,boot_table);
    for b=1:599; bootF(b,:) = model.continuous.F{b}(1,:); end
    p_value1_OLS_boot(MC) = 1 - (sum(F1_obs>bootF(:,1)) / 599);
    p_value2_OLS_boot(MC) = 1 - (sum(F2_obs>bootF(:,2)) / 599);

    % if we use IRLS
    % ---------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'IRLS');
    F1_obs =   model.continuous.F(1,1);
    F2_obs =   model.continuous.F(1,2);
    p_value1_IRLS(MC) = model.continuous.p(1,1);
    p_value2_IRLS(MC) = model.continuous.p(1,2);
    
    % if we use IRLS and bootstrap
    % ----------------------------
    weights = model.W';
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,weights,boot_table);
    for b=1:599; bootF(b,:) = model.continuous.F{b}(1,:); end
    p_value1_IRLS_boot(MC) = 1 - (sum(F1_obs>bootF(:,1)) / 599);
    p_value2_IRLS_boot(MC) = 1 - (sum(F2_obs>bootF(:,2)) / 599);
end
save multiple_regression;    


%% ANOVA
% we compare OLS to IRLS with/without bootstrap
% case 1: data are normally distributed 
clear

% N=3*20 3 groups
% -----------------------
for MC = 1:1000
    r = randn(20,3);
    Y(1,:) = r(:);
    Y(2,:) = r(:);
    X = [kron(eye(3),ones(20,1)) ones(60,1)];
    nb_conditions = 3;
    nb_continuous = 0;
    nb_interactions = 0;
    
    % if we use OLS
    % --------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'OLS');
    F_obs =  model.conditions.F(1);
    p_value_OLS(MC) = model.conditions.p(1);
    
    % if we use OLS and bootstrap
    % --------------------------
    weights = model.W';
    boot_table = randi(60,60,599);
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,weights,boot_table);
    for b=1:599; bootF(b) =  model.conditions.F{b}(1); end
    p_value_OLS_boot(MC) = 1 - (sum(F_obs>bootF) / 599);
    
    % if we use IRLS
    % ---------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'IRLS');
    F_obs =  model.conditions.F(1);
    p_value_IRLS(MC) = model.conditions.p(1);
    
    % if we use IRLS and bootstrap
    % ----------------------------
    weights = model.W';
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,weights,boot_table);
    for b=1:599; bootF(b) =  model.conditions.F{b}(1); end
    p_value_IRLS_boot(MC) = 1 - (sum(F_obs>bootF) / 599);
end
save ANOVA1;    


clear
% N=4*20 2*2 ANOVA with interaction
% --------------------------------
for MC = 1:1000
    r = randn(20,4);
    Y(1,:) = r(:);
    Y(2,:) = r(:);
    X = [kron(eye(2),ones(40,1)) repmat(kron(eye(2),ones(20,1)),2,1)];
    nb_conditions = [2 2];
    nb_continuous = 0;
    [tmpX nb_interactions] = limo_make_interactions(X, nb_conditions);
    X = [tmpX ones(80,1)];
    
    % if we use OLS
    % --------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'OLS');
    F1_obs =  model.conditions.F(1,1);
    F2_obs =  model.conditions.F(2,1);
    F3_obs =  model.interactions.F(1);
    p_value1_OLS(MC) = model.conditions.p(1,1);
    p_value2_OLS(MC) = model.conditions.p(2,1);
    p_value3_OLS(MC) = model.interactions.p(1);

    % if we use OLS and bootstrap
    % --------------------------
    weights = model.W';
    boot_table = randi(80,80,599);
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,weights,boot_table);
    for b=1:599; bootF(b,1:2) =  model.conditions.F{b}(:,1); bootF(b,3) =  model.interactions.F{b}(1); end
    p_value1_OLS_boot(MC) = 1 - (sum(F1_obs>bootF(:,1)) / 599);
    p_value2_OLS_boot(MC) = 1 - (sum(F2_obs>bootF(:,2)) / 599);
    p_value3_OLS_boot(MC) = 1 - (sum(F3_obs>bootF(:,3)) / 599);
    
    
    % if we use IRLS
    % ---------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'IRLS');
    F1_obs =  model.conditions.F(1,1);
    F2_obs =  model.conditions.F(2,1);
    F3_obs =  model.interactions.F(1);
    p_value1_IRLS(MC) = model.conditions.p(1,1);
    p_value2_IRLS(MC) = model.conditions.p(2,1);
    p_value3_IRLS(MC) = model.interactions.p(1);
    
    % if we use IRLS and bootstrap
    % ----------------------------
    weights = model.W';
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,weights,boot_table);
    for b=1:599; bootF(b,1:2) =  model.conditions.F{b}(:,1); bootF(b,3) =  model.interactions.F{b}(1); end
    p_value1_OLS_boot(MC) = 1 - (sum(F1_obs>bootF(:,1)) / 599);
    p_value2_OLS_boot(MC) = 1 - (sum(F2_obs>bootF(:,2)) / 599);
    p_value3_OLS_boot(MC) = 1 - (sum(F3_obs>bootF(:,3)) / 599);
end
save ANOVA2;    


%% ANCOVA
% we compare OLS to IRLS with/without bootstrap
% case 1: data are normally distributed 
clear

% N=3*20 3 groups and 1 cov
% -------------------------
for MC = 1:1000
    r = randn(20,3);
    cov = randn(60,1);
    Y(1,:) = r(:);
    Y(2,:) = r(:);
    X = [kron(eye(3),ones(20,1)) ones(60,1)];
    cov = cov - X*pinv(X)*cov; % cov is now orthogonal to groups
    X = [kron(eye(3),ones(20,1)) cov ones(60,1)];
    nb_conditions = 3;
    nb_continuous = 1;
    nb_interactions = 0;
    
    % if we use OLS
    % --------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'OLS');
    F1_obs =  model.conditions.F(1);
    F2_obs =  model.continuous.F(1);
    p_value1_OLS(MC) = model.conditions.p(1);
    p_value2_OLS(MC) = model.continuous.p(1);
    
    % if we use OLS and bootstrap
    % --------------------------
    weights = model.W';
    boot_table = randi(60,60,599);
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,weights,boot_table);
    for b=1:599; bootF(b,1) =  model.conditions.F{b}(1); bootF(b,2) =  model.continuous.F{b}(1); end
    p_value1_OLS_boot(MC) = 1 - (sum(F1_obs>bootF(:,1)) / 599);
    p_value2_OLS_boot(MC) = 1 - (sum(F2_obs>bootF(:,2)) / 599);

    % if we use IRLS
    % ---------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'IRLS');
    F1_obs =  model.conditions.F(1);
    F2_obs =  model.continuous.F(1);
    p_value1_IRLS(MC) = model.conditions.p(1);
    p_value2_IRLS(MC) = model.continuous.p(1);
    
    % if we use IRLS and bootstrap
    % ----------------------------
    weights = model.W';
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,weights,boot_table);
    for b=1:599; bootF(b) =  model.F{b}(1); end
    for b=1:599; bootF(b,1) =  model.conditions.F{b}(1); bootF(b,2) =  model.continuous.F{b}(1); end
    p_value1_IRLS_boot(MC) = 1 - (sum(F1_obs>bootF(:,1)) / 599);
    p_value2_IRLS_boot(MC) = 1 - (sum(F2_obs>bootF(:,2)) / 599);
end
save ANCOVA1;    

clear

% N=3*20 3 groups and 1 cov
% -------------------------
for MC = 1:1000
    r = randn(20,3) + repmat([0 5 10],20,1);
    cov = randn(60,1); 
    Y(1,:) = r(:);
    Y(2,:) = r(:);
    X = [kron(eye(3),ones(20,1)) ones(60,1)];
    cov = cov - X*pinv(X)*cov; % cov is now orthogonal to groups
    cov(1:20) = cov(1:20) - mean(cov(1:20));
    cov(21:40) = cov(21:40) + 5;
    cov(41:60) = cov(41:60) + 10;
    X = [kron(eye(3),ones(20,1)) cov ones(60,1)];
    nb_conditions = 3;
    nb_continuous = 1;
    nb_interactions = 0;
    
    % if we use OLS
    % --------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'OLS');
    F1_obs =  model.conditions.F(1);
    F2_obs =  model.continuous.F(1);
    p_value1_OLS(MC) = model.conditions.p(1);
    p_value2_OLS(MC) = model.continuous.p(1);
    
    % if we use OLS and bootstrap
    % --------------------------
    weights = model.W';
    boot_table = randi(60,60,599);
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,weights,boot_table);
    for b=1:599; bootF(b,1) =  model.conditions.F{b}(1); bootF(b,2) =  model.continuous.F{b}(1); end
    p_value1_OLS_boot(MC) = 1 - (sum(F1_obs>bootF(:,1)) / 599);
    p_value2_OLS_boot(MC) = 1 - (sum(F2_obs>bootF(:,2)) / 599);

    % if we use IRLS
    % ---------------
    model = limo_glm1(Y',X,nb_conditions,nb_interactions,nb_continuous,'IRLS');
    F1_obs =  model.conditions.F(1);
    F2_obs =  model.continuous.F(1);
    p_value1_IRLS(MC) = model.conditions.p(1);
    p_value2_IRLS(MC) = model.continuous.p(1);
    
    % if we use IRLS and bootstrap
    % ----------------------------
    weights = model.W';
    model = limo_glm1_boot(Y',X,nb_conditions,nb_interactions,nb_continuous,1,weights,boot_table);
    for b=1:599; bootF(b) =  model.F{b}(1); end
    for b=1:599; bootF(b,1) =  model.conditions.F{b}(1); bootF(b,2) =  model.continuous.F{b}(1); end
    p_value1_IRLS_boot(MC) = 1 - (sum(F1_obs>bootF(:,1)) / 599);
    p_value2_IRLS_boot(MC) = 1 - (sum(F2_obs>bootF(:,2)) / 599);
end
save ANCOVA2; 



