function limo_test_AIC_BIC

% unit test function for limo_AIC_BIC.m
% 
% providing known input parameters, the functions must returns known model parameters
% the known values (ground truth) values are obtained using matlab fitglm (stat toolbox)
% 

mkdir('tmp'); cd('tmp'); 
directory = pwd;

% we make data as for EEG (same data used in limo_test_glm.m)
% ----------------------------------------------------------
Y        = NaN(1,1,15);
Y(1,1,:) = [5,6,8,7,9,3,2,1,5,6,4,5,8,9,6];
Cat      = [1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 ; 1 1 2 2 3 1 1 2 2 3 3 3 1 2 3];
Cont     = [0.1978 1.3107 0.5688 -0.5441 -1.286 -0.915 -0.1731 0.1978 0.5688 0.9398 1.6817 -0.915 0.9398 -1.286 -1.286 ; -1.0185 -0.3542 0.31 -0.3542 -1.0185 0.9742 -0.3542 0.31 2.3026 1.6384 -0.3542 -1.0185 -1.0185  -0.3542 0.31]';
[X,nb_conditions,nb_interactions,nb_continuous] = ...
    limo_design_matrix(Y,Cat,Cont,directory,1,0,1);
title('GLM with 9 conditions and 2 covariates (3*3 ANCOVA)')

model1     = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'OLS','Time',[],[]);
[aic, bic] = limo_AIC_BIC('Yr',squeeze(Y(1,1,:)),'Betas',model1.betas,'X',X);
% https://se.mathworks.com/help/stats/fitlm.html
X2         = X;
X2(find(X(:,3)),[1 2]) = -1;
X2(find(X(:,6)),[4 5]) = -1;
X2(:,[3 6 9]) = [];
model2     = fitlm(X2,squeeze(Y(1,1,:)));
if model2.ModelCriterion.BIC == bic(1) && model2.ModelCriterion.AIC == aic(1)
    disp('-------')
    disp('success')
    disp('-------')
else
    error('LIMO AIC and/or BIC diverge from fitlm')
end

cd ..; try rmdir('tmp','s'); end

