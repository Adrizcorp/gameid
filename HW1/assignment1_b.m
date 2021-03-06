%% EECE 574 Self-tuning
%% Author: Holguer A. Becerra
%% assignment 1.
%% Professor: Guy Dummont.
%% RLS Second order a and b

clear
clc
load data2017.mat

u=data(:,1); %System Input
y=data(:,2); %System Output

n=2;% system order
m=length(u);

theta=zeros(1,2*n)'; % create a zero vector for the coefficients a and b, 2 for a 2 for b
alpha=1e4; %%init factor
P=alpha*eye(2*n); % Covariance Matrix
lambda=1; % Forgetting Factor
lambda_inv=1/lambda; % inverse Forgetting Factor

    for k=3:m,%%sweeping out y
        %taking the new seeds for the coefficients estimation or observations,
        %2 samples behind
        xt=[-y(k-1) -y(k-2) u(k-1) u(k-2)];
        x=xt';%this would be xT-- or x Transpose
        %calculating P cvariance matrix
        P=lambda_inv*(P-(P*x*xt*P)/(lambda+xt*P*x));
        %taking the covariance matrix calculate and update the weights (Thetha, which are A(t), and B(t) coeffiecients)
        theta=theta-P*x*(xt*theta-y(k));
                   %____
                     %|
                     %.-> K(t+1)=P*x=P(t+1)*x(t+1)
    end

%taking out the coeficcients
a1=theta(1,1);
a2=theta(2,1);
b1=theta(3,1);
b2=theta(4,1);
u1=u;%copy of the input signal
numerator=[b1 b2];
denomi=[1 a1 a2];
yestimate=dlsim(numerator,denomi,u1);%simulation of a discrete linear system to calculate Y estimated

%%
corre_value= correlate_signals(y,yestimate);%%calculate the correlation between the 2 signals, the higher the value the
%more correlated or similar the signals are.
str=sprintf('similarity= %f%%', corre_value);
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


