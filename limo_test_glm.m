function limo_test_glm

% unit test function for limo_glm with OLS estimation +
% limo_semi_partial_coef and limo_contrast - ensuring all GLM effects are
% well estimated
% 
% providing known input parameters, the functions must returns known model parameters
% the known values (ground truth) values were obtained using JASP
%
%   model = limo_glm(Y,X,nb_conditions,nb_interactions,nb_continuous, method,analysis type,n_freqs,n_times)
%
%   model.R2_univariate = the R2 of the model
%   model.F             = the F value of the model
%   model.df            = the df associated to the model F
%   model.p             = the p value of the model
%   model.betas dim     = the beta parameters (dimensions nb of paramters x frames)
%   model.conditions    = categorical effects
%          --> F/p in rows are the factors, in columns time frames
%          --> df row 1 = df, row 2 = dfe, columns are factors
%   model.continuous = continuous effects
%          --> F/p in rows are the variables, in columns time frames
%          --> df column 1 = df, column2 2 = dfe (same for all covariates)
%
% Cyril R. Pernet
% University of Edinburgh
% -----------------------
index = 1; errorcatch = {};
mkdir('tmp'); cd('tmp'); 
directory = pwd;

%% 0 groups - just the mean 
% (useful for power analyes when using WLS)
Y        = NaN(1,1,10);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6]';
Cat      = [];
Cont     = [];
[X,nb_conditions,nb_interactions,nb_continuous] = ...
    limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
title('GLM with just the mean')
model = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if single(model.betas) == single(mean(Y))
    disp('----------------------------')    
    disp('GLM no regressors validated ')
    disp('----------------------------')    
else
    errorcatch{index} = 'GLM just the mean failed';
    index = index + 1; disp('GLM with just the mean failed')
end


%% 2 groups
Y        = NaN(1,1,10);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6]';
Cat      = [1 1 1 1 1 2 2 2 2 2]';
Cont     = [];
[X,nb_conditions,nb_interactions,nb_continuous] = ...
    limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
title('GLM with two conditions (ANOVA)')
model = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if (model.conditions.F - 9.529412) < 0.0001 && ...
        (model.conditions.df(1) -1) == 0 && ...
        (model.conditions.df(2)-8) == 0 && ...
        (model.conditions.p - 0.014958)< 0.0001
    disp('----------------------------')    
    disp('ANOVA two samples validated')
    disp('----------------------------')    
else
    errorcatch{index} = 'ANOVA two samples failed';
    index = index + 1; disp('ANOVA two samples failed')
end

%% 1 way ANOVA

Y        =  NaN(1,2,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Y(1,2,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat      = [1 1 1 1 1 2 2 2 2 2 2 2 3 3 3];
Cont     = 0;
[X,nb_conditions,nb_interactions,nb_continuous] = ...
    limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
title('GLM with three conditions (ANOVA)')
model = limo_glm(squeeze(Y(1,:,:))',X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
load('Betas.mat'); Betas(1,:,:) = model.betas'; save ('Betas.mat','Betas')
load('Res.mat'); Res(1,:,:)  = (squeeze(Y(1,:,:))' - X*model.betas)'; save('Res.mat','Res')
if (model.R2_univariate(1) - 0.582166) < 0.001 && ...
        (model.conditions.F(1) - 8.359777) < 0.001 && ...
        (model.conditions.df(1) -2) < 0.001 && ...
        (model.conditions.df(2)-12) < 0.001 && ...
        (model.conditions.p(1) - 0.005321) < 0.001
    disp('----------------------')
    disp('1 way ANOVA validated')
    disp('----------------------')    
else
    errorcatch{index} = '1 way ANOVA failed';
    index = index + 1; disp('1 way ANOVA failed')
end

% compute pair-wise contrasts
LIMO.dir                     = pwd;
LIMO.Type                    = 'Channels';
LIMO.Analysis                = 'Time';
LIMO.Level                   = 1;
LIMO.design.X                = X;
LIMO.design.type_of_analysis = 'Mass-univariate';
LIMO.design.method           = 'OLS';
LIMO.model.model_df          = model.df; 
LIMO.contrast{1}.C           = limo_contrast_checking(LIMO.dir, LIMO.design.X, [1 -1]);
LIMO.contrast{1}.V           = 'T';
save LIMO LIMO
result1                      = limo_contrast(Y, Betas, LIMO, 'T', 1);
LIMO.contrast{2}.C           = limo_contrast_checking(LIMO.dir, LIMO.design.X, [0 1 -1]);
LIMO.contrast{2}.V           = 'T';
save LIMO LIMO
result2                      = limo_contrast(Y, Betas, LIMO, 'T', 1);
if (result1(1,1,4) - 3.32902) < 0.001 && ...
        (result1(1,1,5) - 0.006009) < 0.001 && ...
        (result2(1,1,4) + 3.39791) < 0.001 && ...
        (result2(1,1,5) - 0.005290) < 0.001 
    disp('----------------------')    
    disp('contrast T ANOVA validated')
    disp('----------------------')    
else
    errorcatch{index} = 'contrast T ANOVA failed';
    index = index + 1; disp('contrast T ANOVA failed')
end

%% 2 ways ANOVA 
Y        = NaN(1,2,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Y(1,2,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat      = [1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 ; 1 1 2 2 3 1 1 2 2 3 3 3 1 2 3];
Cont     = 0;
[X,nb_conditions,nb_interactions,nb_continuous] = ...
    limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
title('GLM with nine conditions (3 * 3 ANOVA)')
model = limo_glm(squeeze(Y(1,:,:))',X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
load('R2.mat'); R2(1,:,1) = model.R2_univariate;
R2(1,:,2) = model.F; R2(1,:,3) = model.p; save ('R2.mat','R2')
if (model.R2_univariate(1) - 0.6939) < 0.001 && ...
        (model.F(1) - 5.668217) < 0.001 && ...
        sum(model.df -[4 10]) == 0 && ...
        (model.p(1) - 0.012004) < 0.001 && ...
        sum(model.conditions.F(:,1)' - [10.3755 1.8259]) < 0.001 && ...
        sum(single(model.conditions.df(1,:)) - [2 10]) == 0 && ...
        sum(model.conditions.p(:,1)' - [0.003637 0.210885]) < 0.001
    disp('----------------------')
    disp('2 ways ANOVA validated')
    disp('----------------------')
else
    errorcatch{index} = '2 ways ANOVA failed';
    index = index + 1; disp('2 ways ANOVA failed')
end


%% %% 2 ways ANOVA with an interaction

 Y        = NaN(1,2,12);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5];
Y(1,2,:) = [5,6,8,7,9,3,2,1,5,6,4,5];
Cat      = [1 1 1 1 1 1 2 2 2 2 2 2 ; 1 2 3 1 2 3 1 2 3 1 2 3 ];
Cont     = 0;
[X,nb_conditions,nb_interactions,nb_continuous] = ...
    limo_design_matrix(Y,Cat,Cont,directory,1,1,1); % note the flag for interaction is 1
title('GLM with six conditions (2 * 3 ANOVA with interaction)')
load Yr % here because of the factor structure Y is reorganized 
model = limo_glm(squeeze(Yr(1,:,:))',X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
load('Betas.mat'); Betas(1,:,:) = model.betas'; save ('Betas.mat','Betas')
load('Res.mat'); Res(1,:,:)  = (squeeze(Y(1,:,:))' - X*model.betas)'; save('Res.mat','Res')
if (model.R2_univariate(1) - 0.482900) < 0.001 && ...
        (model.F(1) - 1.120635) < 0.001 && ...
        sum(model.df -[5 6]) == 0 && ...
        (model.p(1) -0.438895) < 0.001 && ...
        sum(model.conditions.F(:,1)' - [3.57143 0.01587]) < 0.001 && ...
        sum(single(model.conditions.df(1,:)) - [1 6]) == 0 && ...
        sum(single(model.conditions.df(2,:)) - [2 6]) == 0 && ...
        sum(model.conditions.p(:,1)' - [0.107679 0.984293]) < 0.001 && ...
        single(model.interactions.F(1)) - 1 == 0 && ...
        sum(single(model.interactions.df) - [2 6]) == 0 && ...
        sum(model.interactions.p(:,1) - 0.421875) < 0.001
    disp('---------------------------------------')
    disp('2 ways ANOVA with interaction validated')
    disp('---------------------------------------')
else
    errorcatch{index} = '2 ways ANOVA with interaction failed';
    index = index + 1; disp('2 ways ANOVA with interaction failed')
end

%% 3 ways ANOVA

Y        = NaN(1,1,24);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6,8,7,6,4,2,3,6,5,9];
F1       = [repmat(1,1,8) repmat(2,1,8) repmat(3,1,8)];
F2       = repmat([repmat(1,1,4) repmat(2,1,4)],1,3);
F3       = repmat([1 1 2 2],1,6);
Cat      = [F1;F2;F3];
Cont     = 0;
[X,nb_conditions,nb_interactions,nb_continuous] = ...
    limo_design_matrix(Y,Cat,Cont,directory,1,1,1);
title('GLM with twelve conditions (3 * 2 * 2 ANOVA)')
model = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if (model.R2_univariate(1) - 0.693136) < 0.001 && ...
        (model.F(1) - 2.464115) < 0.001 && ...
        sum(single(model.df) -[11 12]) == 0 && ...
        (model.p(1) -0.068242) < 0.001 && ...
        sum(model.conditions.F' - [1.1974 0.2105 1.8947]) < 0.001 && ...
        sum(single(model.conditions.df(:)) - [2 1 1 12 12 12]') == 0 && ...
        sum(model.conditions.p' - [0.335633 0.654555 0.193810]) < 0.001 && ...
        sum(model.interactions.F' - [4.9868 0.1184 0.0526 6.1711 ]) < 0.001 && ...
        sum(single(model.interactions.df(:,1))' - [2 2 1 2]) == 0 && ...
        sum(model.interactions.p' - [0.026526 0.889347 0.822409 0.014353 ])< 0.001 
    disp('---------------------------------------')
    disp('3 ways ANOVA with interaction validated')
    disp('---------------------------------------')
else
    errorcatch{index} = '3 ways ANOVA with interaction failed';
    index = index + 1; disp('3 ways ANOVA with interaction failed')
end

%% simple regression
Y        = NaN(1,1,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat      = [];
Cont     = [0.1978 1.3107 0.5688 -0.5441 -1.286 -0.915 -0.1731 0.1978 0.5688 0.9398 1.6817 -0.915 0.9398 -1.286 -1.286];
[X,nb_conditions,nb_interactions,nb_continuous] = ...
    limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
title('GLM with 1 continuous variable (simgple regression)')
model = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if (model.R2_univariate - 0.031613) < 0.001 && ... 
    (model.continuous.F - 0.424388)  < 0.001 && ...
        sum(single(model.continuous.df) - [1 13]) == 0 && ...
        (model.continuous.p - 0.526105) < 0.001
    disp('---------------------------')    
    disp('Simple regression validated')
    disp('---------------------------')    
else
    errorcatch{index} = 'Simple regression failed';
    index = index + 1; disp('Simple regression failed')
end

%% multiple regression

Y        = NaN(1,2,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Y(1,2,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat      = [];
Cont     = [0.1978 1.3107 0.5688 -0.5441 -1.286 -0.915 -0.1731 0.1978 0.5688 0.9398 1.6817 -0.915 0.9398 -1.286 -1.286 ; ...
    -1.0185 -0.3542 0.31 -0.3542 -1.0185 0.9742 -0.3542 0.31 2.3026 1.6384 -0.3542 -1.0185 -1.0185  -0.3542 0.31]';
[X,nb_conditions,nb_interactions,nb_continuous] = ...
    limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
title('GLM with 2 continuous variables (multiple regression)')
model = limo_glm(squeeze(Y)',X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
load('R2.mat'); R2(1,:,1) = model.R2_univariate;
R2(1,:,2) = model.F; R2(1,:,3) = model.p; save ('R2.mat','R2')
if (model.F(1) - 0.531554332) < 0.001 && ...
        sum(single(model.df) - [2 12]) == 0 && ...
        (model.p(1) - 0.6009) < 0.001 && ...
        sum(model.continuous.F(1,:) - [0.2363 0.6501])< 0.001 && ...
        sum(model.continuous.p(1,:) - [0.6355 0.4358])< 0.001 && ...
        sum(single(model.continuous.df) - [1 12])==0
    disp('-----------------------------')    
    disp('Multiple regression validated')
    disp('-----------------------------')    
else
    errorcatch{index} = 'Multiple regression failed';
    index = index + 1; disp('Multiple regression failed')
end

% compute semi-partial coef
limo_semi_partial_coef(pwd,'Time','Channels',X,model.W,...
    nb_conditions,nb_interactions,nb_continuous,'OLS',model.df)
SPC1 = load('semi_partial_coef_1.mat');
SPC2 = load('semi_partial_coef_2.mat');
if  (SPC1.semi_partial_coef(1,2) - (-0.48634711)^2) < 0.001 && ...
        (SPC2.semi_partial_coef(1,2) - (-0.806313503)^2) < 0.001  
    disp('----------------------')    
    disp('semi partial coefficient validated')
    disp('----------------------')    
else
    errorcatch{index} = 'semi partial coefficient failed';
    index = index + 1; disp('semi partial coefficient failed')
end

%% 1 way ANCOVA
Y        = NaN(1,1,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat      = [1 1 1 1 1 2 2 2 2 2 2 2 3 3 3];
Cont     = [0.1978 1.3107 0.5688 -0.5441 -1.286 -0.915 -0.1731 0.1978 0.5688 0.9398 1.6817 -0.915 0.9398 -1.286 -1.286 ; -1.0185 -0.3542 0.31 -0.3542 -1.0185 0.9742 -0.3542 0.31 2.3026 1.6384 -0.3542 -1.0185 -1.0185  -0.3542 0.31]';
[X,nb_conditions,nb_interactions,nb_continuous] = ...
    limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
title('GLM with 3 conditions and 2 covariates (ANCOVA)')
model = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if (model.F - 3.7277) < 0.001 && ...
        sum(single(model.df) -[4 10]) == 0 && ...
        (model.p - 0.0416) < 0.001 && ...
        (model.conditions.F - 6.4416) < 0.001 && ...
        (model.conditions.p - 0.0159) < 0.001 && ...
        sum(single(model.conditions.df)' - [2 10]) == 0 && ...
        sum(model.continuous.F - [0.0174 0.4052])< 0.001 && ...
        sum(model.continuous.p - [0.8978 0.5387])< 0.001 && ...
        sum(single(model.continuous.df) - [1 10]) ==  0
    disp('----------------------')    
    disp('1 way ANCOVA validated')
    disp('----------------------')    
else
    errorcatch{index} = '1 way ANCOVA failed';
    index = index + 1; disp('1 way ANCOVA failed')
end

%% 2 ways ANCOVA

Y        = NaN(1,1,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat      = [1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 ; 1 1 2 2 3 1 1 2 2 3 3 3 1 2 3];
Cont     = [0.1978 1.3107 0.5688 -0.5441 -1.286 -0.915 -0.1731 0.1978 0.5688 0.9398 1.6817 -0.915 0.9398 -1.286 -1.286 ; -1.0185 -0.3542 0.31 -0.3542 -1.0185 0.9742 -0.3542 0.31 2.3026 1.6384 -0.3542 -1.0185 -1.0185  -0.3542 0.31]';
[X,nb_conditions,nb_interactions,nb_continuous] = ...
    limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
title('GLM with 9 conditions and 2 covariates (3*3 ANCOVA)')
load Yr; model = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if (model.F - 3.214315) < 0.001 && ...
        sum(single(model.df) -[4 10]) == 0 && ...
        (model.p - 0.065197) < 0.001 && ...
        sum(model.conditions.F' - [7.2751 1.4768]) < 0.001 && ...
        sum(model.conditions.p' - [0.015840 0.284533]) < 0.001 && ...
        sum(single(model.conditions.df(1,:)) - [2 8]) == 0 && ...
        sum(model.continuous.F - [0.0593 0.2401])< 0.001 && ...
        sum(model.continuous.p - [0.813750 0.637294])< 0.001 && ...
        sum(single(model.continuous.df) - [1 8]) ==  0
    disp('-----------------------')    
    disp('2 ways ANCOVA validated')
    disp('-----------------------')    
else
    errorcatch{index} = '2 ways ANCOVA failed';
    index = index + 1; disp('2 ways ANCOVA failed')
end

cd ..; try rmdir('tmp','s'); end
if isempty(errorcatch)
    disp('SUCCESS all tests passed !!'); 
else
    cellfun(@(x) fprintf('%s\n',x), errorcatch);
end
