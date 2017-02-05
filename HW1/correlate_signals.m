
%% EECE 574 Self-tuning
%% Author: Holguer A. Becerra
%% assignment 1.
%% Professor: Guy Dummont.
%% Correlation Value

function [fit,yp] = correlate_signals(signal_1,signal_2)
correlation=signal_1*signal_2;
correlation=correlation./sqrt((signal_1.^2)*(signal_2.^2));

fit=0;
LEN=length(signal_1);
for i = 1: LEN
      fit=signal_1(i)*signal_2(i)+fit;
end

yp=correlation;


