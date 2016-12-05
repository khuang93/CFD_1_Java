clear all;
close all
clc;
M=csvread('plot_data_pr3_t_4000_java.csv',1,0);
% M=csvread('plot_data_pr3_final_java.csv',1,0);

% M=csvread('solData.csv',1,0);
sizeM=size(M);
numRows=sizeM(1);
numCols=sizeM(2);
x=M(1:numRows,1:1);
Ts=M(1:numRows,2:2);
Tf=M(1:numRows,3:3);

plot(x,Ts);
hold on
plot(x,Tf);
legend Ts Tf
xlabel 'x in m'
ylabel 'T in °C'
%%
TT=csvread('plot_cycles.csv',0,0);
sizeTT=size(TT);
numRowsTT=sizeTT(1);
numColsTT=sizeTT(2);
time=TT(1:numRowsTT,1:1);
status=TT(1:numRowsTT,2:2);
plot(time, status);