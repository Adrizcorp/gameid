%Author: Alireza Taale
%Date: FEB 2,2017
%Purpose: System identification HW1
%Input Excitation: given
%UBC Electrical Engineering Department, Vancouver, Canada
%Assignment #1
%File Path:C:\Users\Alireza\Google Drive\CourseWorks\Jan-2017\EECE-574\
%12Homework 1 - Due February 11, 2017
%File Name: ass01_aV2
%Modifications: Preprocessing was added= Tapering by a hann window and zero
%
%************************************************************************%
clear;
clc;
load('data2017');
u = data(:,1); y = data(:,2);
u = detrend(u);
y = detrend(y);
data = iddata(y,u,1,'InputName','u(t)','OutputName','y(t)');
advice(data)
isnlarx(data,[1,2,3]);
ir=cra(data) 
sr=cumsum(ir) 
figure;
plot(sr);

%sys = arx(data, [2 2 2]);
dt = 1;
N = length(u);



%% create an ARX model
estimate = data(1:700);
verif = data(701:1025);
v = arxstruc(estimate,verif,struc(1:10, 1:10, 1:10));
% give me advice about the order I could take
order = selstruc(v,'aic');
% base on that start from 1 to the ARX suggestion order for ARMAX
na = 1:order(1);
nc = 1:order(3);
nk = 0:2;
L = length(na) * length(nc) * length(nk);
models = cell(1,L);
ct = 1;
for i = 1:length(na)
    na_ = na(i);
    nb_ = na_;
    for j = 1:length(nc)
        nc_ = nc(j);
        for k = 1:length(nk)
            nk_ = nk(k); 
            models{ct} = armax(data,[na_ nb_ nc_ nk_]);
            ct = ct+1;
        end
    end
end  

%% 
% Stack the estimated models and compare their simulated responses to estimation
% data |z|. 
models = stack(1,models{:});
compare(data,models)   

%########################################################################%

% Plot the input and output in time domain

%Preprocessing
%Hann filter
u = u .* hann(N);
y = y .* hann(N);
u = [ u ; zeros(2024-N,1)];
y = [ y ; zeros(2024-N,1)];
N = length(y);
%zeropaded

figure();
subplot(2,1,1);
plot(u);
grid on;
ylim([-.2 1.22]);

subplot(2,1,2);
plot(y);
grid on;
ylim([-2 2]);

% plot the input and output in frequency domain

fy = fft(y,N);
fu = fft(detrend(u),N);
faxis= [-N/2:-1 , 0:N/2-1 ] *1 /(N*dt); % for even N
% faxis= [-(N-1)/2:-1 , 0:(N-1)/2]; % for odd N

figure();
subplot(2,2,1);
plot(faxis, fftshift(abs(fu)));
ylabel('|U(f)|');

subplot(2,2,2);
plot(faxis, fftshift(angle(fu)));
ylabel('<U(f)');

subplot(2,2,3);
plot(faxis, fftshift(abs(fy)));
ylabel('|Y(f)|');

subplot(2,2,4);
plot(faxis, fftshift(angle(fy)));
ylabel('<Y(f)');


H = fy.*conj(fu)./(fu.*conj(fu)+0.001);

yEST = real(ifft(H.*fu));
% figure();
% plot(yEST);

% figure();
% [Ryyhat,lags]=xcorr(yEST,y);
% plot(Ryyhat);

figure();
plot(y);
hold on;
plot(yEST);

figure();
plot(faxis, fftshift(abs(H)) );
[pks,locs,w,p] = findpeaks(abs(H)) ;

% %Plot the autocorrelation of the imput and cross correlation of the output
% %and input
% [Ruu,lags,bounds] = autocorr(u,[]);
% figure();
% plot(lags,Ruu);
% 
% [Ryu,lags] = xcorr(u,y);
% plot(lags,Ryu);
% 
% 
% 
% yest = conv(u,h);
% figure();
% plot(yest);
% 
% %Split the data to two parts; one for modeling and the other for
% %verification
% 
% %Option 1: impulse response identification
% 
% %Option 2: Generalized least square
% 
% %Option 3: Classic LSE
% 
% 
% 
% %Compare two models