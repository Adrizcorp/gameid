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

na=3;% system order
nb=4;% system order
nc=600;
m=length(u);

theta=zeros(1,na+nb+nc)'; % create a zero vector for the coefficients a and b, 2 for a 2 for b
alpha=1e4; %%init factor
residual=zeros(1,m);
P=alpha*eye(na+nb+nc); % Covariance Matrix
lambda=1; % Forgetting Factor
lambda_inv=1/lambda; % inverse Forgetting Factor

    for k=(na+nb+nc):m,%%sweeping out y
        %taking the new seeds for the coefficients estimation or observations,
        %2 samples behind
        phit=[];
        for order=1:na,
            phit=[phit -y(k-order)];
        end
        for order=1:nb,
            phit=[phit u(k-order)];
        end
        for order=1:nc,
            phit=[phit residual(k-order)];
        end
        residual(k)=y(k)-phit*theta;
        phi=phit';%this would be xT-- or x Transpose
        %calculating P cvariance matrix
        P=lambda_inv*(P-(P*phi*phit*P)/(lambda+phit*P*phi));
        %taking the covariance matrix calculate and update the weights (Thetha, which are A(t), and B(t) coeffiecients)
        theta=theta-P*phi*(phit*theta-y(k));
                   %____
                     %|
                     %.-> K(t+1)=P*phi=P(t+1)*x(t+1)
    end

%taking out the coeficcients
a=[];
b=[];
c=[];
for order=1:na,
   a=[a theta(order,1)];
end
for order=1:nb,
   b=[b theta(na+order,1)];
end
for order=1:nc,
   c=[c theta(na+nb+order,1)];
end
u1=u;%copy of the input signal
numerator=b;
denomi=[1 a];
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

A=[1 a];
B=[b];
C = [1 c];
arx1red=idpoly(A,B,C,1); 
zi=iddata(y(100+21:502),u(20+100:501),2); 
zi2=iddata(y(500+21:1025),u(20+500:1024),2); 
compare(zi,arx1red);
compare(zi2,arx1red);

