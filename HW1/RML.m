%% EECE 574 Self-tuning
%% Author: Holguer A. Becerra
%% assignment 1.
%% Professor: Guy Dummont.
%% RML

clear
clc
load data2017.mat

u=data(:,1); %System Input
y=data(:,2); %System Output

n=2;% system order
m=length(u);

theta=zeros(1,3*n)'; % create a zero vector for the coefficients a and b, 2 for a 2 for b
theta_prev=theta;
error=zeros(1,m);
residual=zeros(1,m);
alpha=1e3; %%init factor
P=alpha*eye(3*n); % Covariance Matrix
lambda=1; % Forgetting Factor
lambda_inv=1/lambda; % inverse Forgetting Factor
for k=3:m,%%sweeping out y
    %taking the new seeds for the coefficients estimation or observations,
    %2 samples behind
    phit=[-y(k-1) -y(k-2) u(k-1) u(k-2) residual(k-1) residual(k-2)];
    phi=phit';%this would be xT-- or x Transpose
    %error(k)=y(k)-phit*theta_prev;
    residual(k)=y(k)-phit*theta;
    %calculating P cvariance matrix
    P=lambda_inv*(P-(P*phi*phit*P)/(lambda+phit*P*phi));
    %taking the covariance matrix calculate and update the weights (Thetha, which are A(t), and B(t) coeffiecients)
    theta=theta-P*phi*(phit*theta-y(k));
               %____
                 %|
                 %.-> K(t+1)=P*phi=P(t+1)*x(t+1)
end

%taking out the coeficcients
a1=theta(1,1);
a2=theta(2,1);
b1=theta(3,1);
b2=theta(4,1);
r1=theta(5,1);
r2=theta(6,1);
u1=u;%copy of the input signal
numerator=[b1 b2];
denomi=[1 a1 a2];
A=[1 a1 a2];
B=[b1 b2];
C = [1 r1 r2];
arx1red=idpoly(denomi,numerator,[1 r1 r2],1); 
zi=iddata(y(2:501),u(1:500),1); 
compare(zi,arx1red)
yestimate=dlsim(numerator,denomi,u1);%simulation of a discrete linear system to calculate Y estimated

%%
corre_value= correlate_signals(y',yestimate);%%calculate the correlation between the 2 signals, the higher the value the
%more correlated or similar the signals are.
str=sprintf('the correlation value of the signal is= %f', corre_value);
%%plot the results.
figure;
plot(yestimate,'r');
hold on;
plot(y,'b');
grid on;
xlabel('samples number');
ylabel('System Response');
dim = [0.2 0.6 0.3 0.3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
legend('Identified Model','Orginal Model');

