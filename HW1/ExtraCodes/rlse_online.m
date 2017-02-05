function [x,K,P] = rlse_online(aT_k1,b_k1,x,P)
K = P*aT_k1'/(aT_k1*P*aT_k1'+1); %Eq.(2.1.17)
x = x +K*(b_k1-aT_k1*x); %Eq.(2.1.16)
P = P-K*aT_k1*P; %Eq.(2.1.18)