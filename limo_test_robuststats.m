function limo_test_robuststats

%% routine to tests robusts statistical tests 
% data and expected results are taken from Wilcox - Introduction to Robust 
% Estimation and Hypothesis testing ; the various limo functions used 
% implement the procedures described in this book, these are called within
% limo_random_robust.m for group level analyses
%
% this function is called against the code base making sure low level stats
% are correct

error_index = 1;

%% 1 sample t-test - limo_trimci.m 

LSAT_scores = [545 555 558 572 575 576 578 580 594 605 635 651 653 661 666];
% make this a 3D matrix
Data = NaN(1,1,length(LSAT_scores));
Data(1,1,:) = LSAT_scores;
[t,tmdata,trimci,se,~,~,df]=limo_trimci(Data);
if strcmp(sprintf('%1.1f',t),'40.0') && ...
        strcmp(sprintf('%1.1f',tmdata),'596.2') && ...
        strcmp(sprintf('%1.2f',se),'14.92') && ...
        strcmp(sprintf('%1.1f',trimci(1)),'561.8') && ...
        strcmp(sprintf('%1.1f',trimci(2)),'630.6')  && df == 8
    disp('1 sample t-test: limo_trimci validated')
else
    whicherror{error_index} = 'limo_trimci problem';
    error_index = error_index+1;
end
clear LSAT_scores Data t tmdata trimci se df

%% two samples t-test: limo_yuen_ttest.m

ozone = [41 38.4 24.4 25.9 21.9 18.3 13.1 27.3 28.5 -16.9 26 17.4 21.8 15.4 27.4 19.2 22.4 17.7 26 29.4 21.4 26.6 22.7; ...
    10.1 6.1 20.4 7.3 14.3 15.5 -9.9 6.8 28.2 17.9 -9 -12.9 14 6.6 12.1 15.7 39.9 -15.9 56.6 -14.7 44.1 -9 NaN];
% make this as two 3D matrices
Data1 = NaN(1,1,length(ozone));
Data1(1,1,:) = ozone(1,:);
Data2 = NaN(1,1,length(ozone)-1); % note here I make data not the same size 
Data2(1,1,:) = ozone(2,1:end-1);
[Ty,~,~,CI,p]=limo_yuen_ttest(Data1,Data2);
if strcmp(sprintf('%1.1f',Ty),'3.4') && ...
        strcmp(sprintf('%1.4f',p),'0.0037') && ...
        strcmp(sprintf('%1.1f',CI(1)),'5.3') && ...
        strcmp(sprintf('%1.2f',CI(2)),'22.85')
    disp('2 samples t-test: limo_yuen_ttest validated')
else
    whicherror{error_index} = 'limo_yuen_ttest problem';
    error_index = error_index+1;
end
clear ozone Data1 Data2 Ty CI p

% robust 1-way ANOVA: limo_robust_1way_anova.m
gp1     = [1 2 3 4 5 6 7 8 9 10];
gp2     = [2 3 4 5 6 7 8 9 10 11];
gp3     = [5 6 7 8 9 10 11 12 13 14];
Y       = NaN(1,30);
Y(1,:)  = [gp1 gp2 gp3]';
X       = kron(eye(3),ones(10,1));
[F,p,~] = limo_robust_1way_anova(Y,X);
if strcmp(sprintf('%1.2f',F),'2.87') && ...
        strcmp(sprintf('%1.1f',p),'0.1')
    disp('Group ANOVA: limo_robust_1way_anova validated')
else
    whicherror{error_index} = 'limo_robust_1way_anova problem';
    error_index = error_index+1;
end
clear gp1 gp2 gp3 Y X F p

% paired t-test: limo_yuend_ttest.m
cholesterol = [190 201 300 240 280 170 280 250 240 220; ...
    210 210 340 190 260 180 200 220 230 200];
% make this as two 3D matrices
Data1 = NaN(1,1,length(cholesterol));
Data1(1,1,:) = cholesterol(1,:);
Data2 = NaN(1,1,length(cholesterol)); 
Data2(1,1,:) = cholesterol(2,:);
[~,diff,se,CI,p,~,df]=limo_yuend_ttest(Data1,Data2);
if strcmp(sprintf('%1.1f',diff),'26.8') && ...
        strcmp(sprintf('%1.4f',p),'0.1536') && ...
        strcmp(sprintf('%1.2f',se),'15.96') && ...
        strcmp(sprintf('%1.1f',CI(1)),'-14.2') && ...
        strcmp(sprintf('%1.2f',CI(2)),'67.87') && df == 5
    disp('paired t-test: limo_yuend_ttest validated')
else
    whicherror{error_index} = 'limo_yuen_ttest problem';
    error_index = error_index+1;
end
clear cholesterol Data1 Data2 diff se CI p df

% robust regression -- handled for now via limo_glm using IRLS -- there is
% a limo_lowess.m hidden function ; we should really get this done

% robust N-ways ANOVA / N-ways ANCOVA - handled for via limo_glm using IRLS

% robust repeated measure ANOVA - limo_robust_rep_anova.m still working on
% it 

if ~exist('whicherror','var')
    disp('-------------------------------')
    disp('all robust stat tests validated')
    disp('-------------------------------')
else
   fprintf('%g error(s) found \n',length(whicherror))
   for e=1:length(whicherror)
       fprintf('error %g: %s\n',e,whicherror{e});
   end
end
