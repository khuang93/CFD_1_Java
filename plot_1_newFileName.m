clear all;
close all
clc;


M=csvread('data_t_100.0h_tCharge_1.csv',1,0);
% M=csvread('plot_data_pr3_final_java.csv',1,0);

% SOL=csvread('data_t_1000.0h_tCharge_4.csv',1,0);
sizeM=size(M);

numRows=sizeM(1);
numCols=sizeM(2);
x=M(1:numRows,2:2);
Tf=M(1:numRows,3:3);
Ts=M(1:numRows,4:4);
Tf_sol=M(1:numRows,5:5);
Ts_sol=M(1:numRows,6:6);

E_Tf=mean(Tf_-Tf_sol)
E_Ts=mean(Ts-Ts_sol)


plot(x,Tf);
hold on
plot(x,Ts);
hold on
plot(x,Tf_sol);
hold on
plot(x,Tf_sol);
% plot(x,Tf_star);
% hold on
% plot(x,Ts_star);
% hold on
legend Tf Ts Tf_sol Ts_sol
xlabel 'x in m'
ylabel 'T in °C'

% sizeSOL=size(SOL);
% numRowsSOL=sizeSOL(1);
% numColsSOL=sizeSOL(2);
% xSOL=SOL(1:numRowsSOL,2:2);
% TfSOL=SOL(1:numRowsSOL,3:3);
% TsSOL=SOL(1:numRowsSOL,4:4);
% 
% plot(xSOL,TfSOL);
% hold on
% plot(xSOL,TsSOL);
% legend Tf Ts Tf-2500 Ts-2500
% xlabel 'x in m'
% ylabel 'T in °C'
%%
% TT=csvread('plot_cycles.csv',0,0);
% sizeTT=size(TT);
% numRowsTT=sizeTT(1);
% numColsTT=sizeTT(2);
% time=TT(1:numRowsTT,1:1);
% status=TT(1:numRowsTT,2:2);
% plot(time, status);