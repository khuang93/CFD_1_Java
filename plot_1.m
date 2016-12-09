clear all;
close all
clc;


M=csvread('plot_data_pr3_t_15000.0_java.csv',1,0);
% M=csvread('plot_data_pr3_final_java.csv',1,0);

SOL=csvread('plot_data_pr3_t_4500.0_java.csv',1,0);
sizeM=size(M);

numRows=sizeM(1);
numCols=sizeM(2);
x=M(1:numRows,1:1);
Tf=M(1:numRows,2:2);
Ts=M(1:numRows,3:3);
Tf_star=M(1:numRows,4:4);
Ts_star=M(1:numRows,5:5);


plot(x,Tf);
hold on
plot(x,Ts);
hold on
% plot(x,Tf_star);
% hold on
% plot(x,Ts_star);
% hold on
legend Tf Ts
xlabel 'x in m'
ylabel 'T in °C'

sizeSOL=size(SOL);
numRowsSOL=sizeSOL(1);
numColsSOL=sizeSOL(2);
xSOL=SOL(1:numRowsSOL,1:1);
TfSOL=SOL(1:numRowsSOL,2:2);
TsSOL=SOL(1:numRowsSOL,3:3);

plot(xSOL,TfSOL);
hold on
plot(xSOL,TsSOL);
legend Tf Ts Tf-2500 Ts-2500
xlabel 'x in m'
ylabel 'T in °C'
%%
% TT=csvread('plot_cycles.csv',0,0);
% sizeTT=size(TT);
% numRowsTT=sizeTT(1);
% numColsTT=sizeTT(2);
% time=TT(1:numRowsTT,1:1);
% status=TT(1:numRowsTT,2:2);
% plot(time, status);