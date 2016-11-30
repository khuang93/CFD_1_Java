clear all;
close all
clc;
M=csvread('plot_data_pr3_t5000java.csv',1,0);
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