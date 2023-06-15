

% EAD

clear


d1 = reshape(load('../Baseline_EAD_DAD_count.csv'), 7 ,[])';
% d1 = reshape(load('../Baseline+75%kmf_EAD_DAD_count.csv'), 7 ,[])';
d1 = reshape(load('../../OLD/Baseline+50%kmf_EAD_DAD_count.csv'), 7 ,[])';
% d1 = reshape(load('../Baseline+25%kmf_EAD_DAD_count.csv'), 7 ,[])';
% d1 = reshape(load('../ISO_EAD_DAD_count.csv'), 7 ,[])';
% d1 = reshape(load('../ISO+75%kmf_EAD_DAD_count.csv'), 7 ,[])';
% d1 = reshape(load('../ISO+50%kmf_EAD_DAD_count.csv'), 7 ,[])';
% d1 = reshape(load('../ISO+25%kmf_EAD_DAD_count.csv'), 7 ,[])';
x= (d1(:,6)>=1);

num1 = sum(x);


x= (d1(:,7)>=1);

num2 = sum(x);

fprintf('\n%d %d\n\n',num1, num2);



% d2 = reshape(load('Baseline_EAD_count.csv'), 7 ,[])';
% 
% x2= (d2(:,6)>4);
% 
% num = sum(x2)
% 
% % d1 = reshape(load('../Baseline_DAD_count.csv'), 6 ,[])';
% % 
% % x= (d1(:,6)>4);
% % 
% % num = sum(x)
% 
% 
% a = d1(:,6);
% b=d2(:,6);
% 
% nn = (x~=x2)