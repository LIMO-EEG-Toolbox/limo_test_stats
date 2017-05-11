% simulation for WLS

mkdir('tmp'); cd('tmp'); 
directory = pwd;
Y = NaN(1,1,100);Cont = [];
Cat = [ones(50,1) ;ones(50,1).*2];

for n=1:1000
R = randn(100,1); Y(1,1,:) = R-mean(R); % H0 is true
[X,nb_conditions,nb_interactions,nb_continuous] =limo_design_matrix(Y,Cat,Cont,directory,1,0,0);
model = limo_glm(squeeze(Y),X,nb_conditions,nb_interactions,nb_continuous,'IRLS','Time',[],[]);
test(n) = (model.p < 0.05); fprintf('round %g done\n',n);
end

cd ..
delete('tmp')
