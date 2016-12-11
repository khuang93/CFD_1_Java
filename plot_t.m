clear;
clc;
% Res=csvread_all();
for i = 72:-1:1
%     Tf(i)=Res(i).matrix(1500,4);
    M=csvread(['plot_data_pr3_t_' num2str(3600*i) '.0_java.csv'],1,0);
    Tf(i)=M(150,4);
end