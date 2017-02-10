%% EECE 574 Self-tuning
%% Author: Holguer A. Becerra
%% assignment 1.
%% Professor: Guy Dummont.
%% RLS

clear
clc
load data2017.mat

u=data(:,1); %System Input
y=data(:,2); %System Output

n=4;% system order
m=length(u);
error=0;
error_n=1;

theta=zeros(1,2*n)'; % create a zero vector for the coefficients a and b, 2 for a 2 for b
alpha=1e4; %%init factor
P=alpha*eye(2*n); % Covariance Matrix
lambda=1; % Forgetting Factor
lambda_inv=1/lambda; % inverse Forgetting Factor

    for k=n*2:m,%%sweeping out y
        %taking the new seeds for the coefficients estimation or observations,
        %2 samples behind
        xt=[];
        for order=1:n,
            xt=[xt -y(k-order)];
        end
        for order=1:n,
            xt=[xt u(k-order)];
        end
        x=xt';%this would be xT-- or x Transpose
        %calculating P cvariance matrix
        P=lambda_inv*(P-(P*x*xt*P)/(lambda+xt*P*x));
        %taking the covariance matrix calculate and update the weights (Thetha, which are A(t), and B(t) coeffiecients)
        theta=theta-P*x*(xt*theta-y(k));
                   %____
                     %|
                     %.-> K(t+1)=P*x=P(t+1)*x(t+1)
         residual=y(k)-xt*theta;
         error=(error+residual*residual);
         error=error*0.5;
         error_v(error_n)=residual;
         error_n=error_n+1;
    end
figure;
plot(error_v);
%taking out the coeficcients
a=[];
b=[];
for order=1:n,
   a=[a theta(order,1)];
   b=[b theta(n+order,1)];
end
u1=u;%copy of the input signal
numerator=b;
denomi=[1 a];
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

A=[1 a];
B=[b];
%C = [1 r1 r2];
arx1red=idpoly(denomi,numerator);%,[1 r1 r2],1); 
zi=iddata(y(100+21:502),u(20+100:501),2); 
compare(zi,arx1red);
hold on;
