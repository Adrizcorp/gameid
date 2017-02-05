%% EECE 574 Self-tuning
%% Author: Holguer A. Becerra
%% assignment 1.
%% Professor: Guy Dummont.

clear
clc
load data2017.mat

u=data(:,1);
y=data(:,2);

n=2;
m=length(u);

theta=zeros(1,2*n)';
alpha=1e3;
P=alpha*eye(2*n);
lambda=1;
lambda_inv=1/lambda;
for k=3:m,%%sweeping out y
    %taking the new seeds for the coefficients estimation or observations,
    %2 samples behind
    phit=[-y(k-1) -y(k-2) u(k-1) u(k-2)];
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
a1=theta(1,1);
a2=theta(2,1);
b1=theta(3,1);
b2=theta(4,1);
u1=u;%copy of the input signal
numerator=[b1 b2];
denomi=[1 a1 a2];
yestimate=dlsim(numerator,denomi,u1);%simulation of a discrete linear system
figure;
plot(yestimate,'r');
hold on;
plot(y,'b');
grid on;