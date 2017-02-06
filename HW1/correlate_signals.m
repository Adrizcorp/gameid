
%% EECE 574 Self-tuning
%% Author: Holguer A. Becerra
%% assignment 1.
%% Professor: Guy Dummont.
%% Correlation Value

function yp= correlate_signals(signal_1,signal_2)
%% normalization
signal_1=detrend(signal_1);
signal_1_abs=abs(signal_1);
signal_1=signal_1/max(signal_1_abs);
signal_2=detrend(signal_2);
signal_2_abs=abs(signal_2);
signal_2=signal_2/max(signal_2_abs);
%% 

correlation=signal_1.*signal_2;
correlation=correlation./(abs(signal_1)).*(abs(signal_2));
yp=1-abs(sum(correlation)/length(correlation));
yp=yp*100;
% negatives=0;
% for i=1:length(correlation),
%     if(correlation(i)<=0)
%     else
%         negatives=negatives+1;
%     end
% end
% 
% yp=(length(correlation)-negatives)/length(correlation);
% 
% 


