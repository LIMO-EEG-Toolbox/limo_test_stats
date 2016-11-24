function limo_test_glm

% unit test function for limo_glm1
% model = limo_glm(Y,X,nb_conditions,nb_interactions,nb_continuous, method,analysis type,n_freqs,n_times)
% 
% provinding know input parameters, the function must returns known model parameters
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
% the parameters and values were obtained using statistica(r)
%
% Cyril R. Pernet
% University of Edinburgh
% -----------------------
index = 1; errorcatch = {};
mkdir('tmp'); cd('tmp'); 

%% 2 groups
directory = pwd;
Y = NaN(1,1,10);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6]';
Cat = [1 1 1 1 1 2 2 2 2 2]';
Cont = [];
[X,nb_conditions,nb_interactions,nb_continuous] =limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
model = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if round(model.conditions.F - 9.5294) == 0 && ...
        round(model.conditions.df(1) -1) == 0 && ...
        round(model.conditions.df(2)-8) == 0 && ...
        round(model.conditions.p - 0.0150) == 0
    disp('two samples t-test validated')
else
    errorcatch{index} = 'two samples t-test failed';
    index = index + 1; disp('two samples t-test failed')
end

%% ANOVA
Y=NaN(1,1,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat = [1 1 1 1 1 2 2 2 2 2 2 2 3 3 3];
Cont = 0;
[X,nb_conditions,nb_interactions,nb_continuous] =limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
model = limo_glm1(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if round(model.conditions.F - 8.3598) == 0 && ...
        round(model.conditions.df(1) -2) == 0 && ...
        round(model.conditions.df(2)-12) == 0 && ...
        round(model.conditions.p - 0.0053) == 0
    disp('1 way ANOVA validated')
else
    errorcatch{index} = '1 way ANOVA failed';
    index = index + 1; disp('1 way ANOVA failed')
end


% 2 ways
Y=NaN(1,1,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat = [1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 ; 1 1 2 2 3 1 1 2 2 3 3 3 1 2 3];
Cont = 0;
[X,nb_conditions,nb_interactions,nb_continuous] =limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
model = limo_glm1(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if sum(round(model.conditions.F' - [10.3755 1.8256])) == 0 && ...
        sum(round(model.conditions.df(1,:) - [2 10])) == 0 && ...
        sum(round(model.conditions.p' - [0.0036 0.2109])) == 0
    disp('2 ways ANOVA validated')
else
    errorcatch{index} = '2 ways ANOVA failed';
    index = index + 1; disp('2 ways ANOVA failed')
end

% % add the intercation
Y=NaN(1,1,12);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5];
Cat = [1 1 1 1 1 1 2 2 2 2 2 2 ; 1 2 3 1 2 3 1 2 3 1 2 3 ];
Cont = 0;
[X,nb_conditions,nb_interactions,nb_continuous] =limo_design_matrix(Y,Cat,Cont,directory,1,1,1); % note the flag for interaction is 1
load Yr % here because of the factor structure Y is reorganized 
model = limo_glm(squeeze(Yr),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
%       SS	    Df	MS	     F	    p
% V1	18.75	1	18.75	3.57143	0.107679
% V2	0.1667	2	0.0833	0.01587	0.9842
% V1*V2	10.5	2	5.25	1	    0.421875
% Error	31.5	6	5.25

if sum(round(model.conditions.F' - [3.57 0.1])) == 0 && ...
        sum(round(model.conditions.df(1,:) - [1 6])) == 0 && ...
        sum(round(model.conditions.df(2,:) - [2 6])) == 0 && ...
        sum(round(model.conditions.p' - [0.1 0.98])) == 0 && ...
        round(model.interactions.F - 1) == 0 && ...
        sum(round(model.interactions.df - [2 6])) == 0 && ...
        round(model.interactions.p - 0.4219) == 0 
    disp('---------------------------------------')
    disp('2 ways ANOVA with interaction validated')
    disp('---------------------------------------')
else
    errorcatch{index} = '2 ways ANOVA with interaction failed';
    index = index + 1; disp('2 ways ANOVA with interaction failed')
end

%% simple regression
Y=NaN(1,1,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat = [];
Cont = [0.1978 1.3107 0.5688 -0.5441 -1.286 -0.915 -0.1731 0.1978 0.5688 0.9398 1.6817 -0.915 0.9398 -1.286 -1.286];
[X,nb_conditions,nb_interactions,nb_continuous] =limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
model = limo_glm1(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if round(model.continuous.F - 0.4244) == 0 && ...
        round(model.continuous.df(1) -1) == 0 && ...
        round(model.continuous.df(2)-13) == 0 && ...
        round(model.continuous.p - 0.5261) == 0
    disp('---------------------------')    
    disp('Simple regression validated')
    disp('---------------------------')    
else
    errorcatch{index} = 'Simple regression failed';
    index = index + 1; disp('Simple regression failed')
end

%% multiple regression

Y=NaN(1,1,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat = [];
Cont = [0.1978 1.3107 0.5688 -0.5441 -1.286 -0.915 -0.1731 0.1978 0.5688 0.9398 1.6817 -0.915 0.9398 -1.286 -1.286 ; -1.0185 -0.3542 0.31 -0.3542 -1.0185 0.9742 -0.3542 0.31 2.3026 1.6384 -0.3542 -1.0185 -1.0185  -0.3542 0.31]';
[X,nb_conditions,nb_interactions,nb_continuous] =limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
model = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if round(model.F - 0.5315) == 0 && ...
        round(model.df(1) -2) == 0 && ...
        round(model.df(2)-12) == 0 && ...
        round(model.p - 0.6009) == 0 && ...
        sum(round(model.continuous.F - [0.2363 0.6501]))==0 && ...
        sum(round(model.continuous.p - [ 0.6355 0.4358]))==0 && ...
        sum(round(model.continuous.df - [1 12]))==0
    disp('-----------------------------')    
    disp('Multiple regression validated')
    disp('-----------------------------')    
else
    errorcatch{index} = 'Multiple regression failed';
    index = index + 1; disp('Multiple regression failed')
end

%% ANCOVA
Y=NaN(1,1,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat = [1 1 1 1 1 2 2 2 2 2 2 2 3 3 3];
Cont = [0.1978 1.3107 0.5688 -0.5441 -1.286 -0.915 -0.1731 0.1978 0.5688 0.9398 1.6817 -0.915 0.9398 -1.286 -1.286 ; -1.0185 -0.3542 0.31 -0.3542 -1.0185 0.9742 -0.3542 0.31 2.3026 1.6384 -0.3542 -1.0185 -1.0185  -0.3542 0.31]';
[X,nb_conditions,nb_interactions,nb_continuous] =limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
model = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if round(model.F - 3.7277) == 0 && ...
        round(model.df(1) -4) == 0 && ...
        round(model.df(2)-10) == 0 && ...
        round(model.p - 0.0416) == 0 && ...
        round(model.conditions.F - 6.4416)==0 && ...
        round(model.conditions.p - 0.0159)==0 && ...
        sum(round(model.conditions.df' - [2 10]))==0 && ...
        sum(round(model.continuous.F - [0.0174 0.4052]))==0 && ...
        sum(round(model.continuous.p - [0.8978 0.5387]))==0 && ...
        sum(round(model.continuous.df - [1 10]))==0
    disp('----------------------')    
    disp('1 way ANCOVA validated')
    disp('----------------------')    
else
    errorcatch{index} = '1 way ANCOVA failed';
    index = index + 1; disp('1 way ANCOVA failed')
end


% 2 ways
Y=NaN(1,1,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat = [1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 ; 1 1 2 2 3 1 1 2 2 3 3 3 1 2 3];
Cont = [0.1978 1.3107 0.5688 -0.5441 -1.286 -0.915 -0.1731 0.1978 0.5688 0.9398 1.6817 -0.915 0.9398 -1.286 -1.286 ; -1.0185 -0.3542 0.31 -0.3542 -1.0185 0.9742 -0.3542 0.31 2.3026 1.6384 -0.3542 -1.0185 -1.0185  -0.3542 0.31]';
[X,nb_conditions,nb_interactions,nb_continuous] =limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
load Yr
model = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
if sum(round(model.conditions.F' - [7.2751 1.4768])) == 0 && ...
        sum(round(model.conditions.df(1,:) - [2 8])) == 0 && ...
        sum(round(model.conditions.p' - [0.0158 0.2845])) == 0 && ...
        sum(round(model.continuous.F - [0.0593 0.2401])) == 0 && ...
        sum(round(model.continuous.p - [0.8138 0.6373])) == 0 && ...
        sum(round(model.continuous.df - [1 8])) == 0
    disp('-----------------------')    
    disp('2 ways ANCOVA validated')
    disp('-----------------------')    
else
    errorcatch{index} = '2 ways ANCOVA failed';
    index = index + 1; disp('2 ways ANCOVA failed')
end


