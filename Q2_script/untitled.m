clear all; close all; clc;

%% Wing
lambda_wing = 1.14/2.68;
ln_wing = 1.82 + 2.77;
c_root = 2.68;
b_half = 4.19;

S = b_half * c_root * (1+lambda_wing);
A_wing = (b_half * 2)^2 / S;

mean_c = ((2*c_root)/3) * (1 + lambda_wing + lambda_wing^2)/(1+ lambda_wing);

%%Tail
lt = (5.88 + 1.33) - (5.67 + 0.78);
lambda_tail = 0.8/(0.67+0.4);

%%NP
dCl_da_wing = (2*pi)/(1 + (2/A_wing));
dCl_da_tail = (2*pi)/(1 + (2/A_wing));
%h_fromwing =