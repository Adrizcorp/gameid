%% EECE 574 Self-tuning
%% Author: Holguer A. Becerra
%% assignment 1.
%% Professor: Guy Dummont.
%% RLS Configurable

clear
clc
load data2017.mat

u=data(:,1); %System Input
y=data(:,2); %System Output


%% 
na=3;% system order for Ay
nb=2;% system order for Bx
nk=1;% delay
%% size of the SISO
m=length(u);
%% Initialize variables.
theta=zeros(1,na+nb)'; % create a zero vector for the coefficients a and b, 2 for a 2 for b
%% transfer function estimations to optimize
init_a=[-0.5784 0.002871 -0.01826];
init_b=[1.5 -1.475];
theta(1)=init_a(1);
theta(2)=init_a(2);
theta(3)=init_a(3);
theta(4)=init_b(1);
theta(5)=init_b(2);

alpha=1e4; %%init factor
P=alpha*eye(na+nb); % Covariance Matrix
lambda=1; % Forgetting Factor
lambda_inv=1/lambda; % inverse Forgetting Factor
%% for errors and parameters graphs
coeffs_historial=zeros(na+nb,m);
coeffs_historial_error=zeros(na+nb,m);
%% algorithm
for k=(na+nb):m,%%sweeping out y
    %taking the new seeds for the coefficients estimation or observations,
    %2 samples behind
    xt=[];
    for order=1:na,
        xt=[xt -y(k-order)];
    end
    for order=1:nb,
        xt=[xt u(k-order)];
    end
    %% RLS calculations
    x=xt';%x from xT
    %calculating P cvariance matrix
    P=lambda_inv*(P-(P*x*xt*P)/(lambda+xt*P*x));
    %taking the covariance matrix calculate and update the weights (Thetha, which are A(t), and B(t) coeffiecients)
    theta_prev=theta;
    theta=theta-P*x*(xt*theta-y(k));
               %____
                 %|
                 %.-> K(t+1)=P*x=P(t+1)*x(t+1)
    %% save the progress of each coefficients to see the convergence, and
    %error. 
    for order=1:na,
        coeffs_historial(order,k)=theta(order);
        coeffs_historial_error(order,k)=abs(coeffs_historial(order,k)-theta_prev(order));
    end
    for order=1:nb,
        coeffs_historial(order+na,k)=theta(order+na);
        coeffs_historial_error(order+na,k)=abs(coeffs_historial(order+na,k)-theta_prev(order+na));
    end

    %% taking out the coeficcients
    a=[];
    b=[];
    for order=1:na,
       a=[a theta(order,1)];
    end
    for order=1:nb,
       b=[b theta(na+order,1)];
    end

end

%% Graph Parameters
texto=[];
figure;
subplot(2,2,[1,2])
for order=1:(na),
    plot(coeffs_historial(order,:));
    texto= [texto;sprintf('a%d',order)];
    hold on;
end
legend(texto);
xlabel('k')
ylabel('parameter')
texto=[];
subplot(2,2,[3,4])
for order=1:(nb),
    plot(coeffs_historial(order+na,:));
    texto= [texto;sprintf('b%d',order)];
    hold on;
end
legend(texto);
xlabel('k')
ylabel('parameter')


%% GRAPH Errors
texto=[];
figure;
subplot(2,2,[1,2])
for order=1:(na),
    plot(coeffs_historial_error(order,:));
    texto= [texto;sprintf('a%d',order)];
    pointText=['min a' num2str(order) '=' num2str(coeffs_historial_error(order,m))];
    %text(m,coeffs_historial_error(order,m)+0.1*order,pointText,'HorizontalAlignment','right');
    hold on;
end

legend(texto);
xlabel('k')
ylabel('Error')
texto=[];
subplot(2,2,[3,4])
for order=1:(nb),
    plot(coeffs_historial_error(order+na,:));
    texto= [texto;sprintf('b%d',order)];
    pointText=['min b' num2str(order) '=' num2str(coeffs_historial_error(order+na,m))];
    %text(m,coeffs_historial_error(order+na,m)+0.1*order,pointText,'HorizontalAlignment','right');
    hold on;
end
legend(texto);
xlabel('k')
ylabel('Error')

%% taking out the coeficcients
a=[];
b=[];
for order=1:na,
   a=[a theta(order,1)];
end
for order=1:nb,
   b=[b theta(na+order,1)];
end
   
u1=u;%copy of the input signal
numerator=b;
denomi=[1 a];
yestimate=dlsim(numerator,denomi,u1);%simulation of a discrete linear system to calculate Y estimated

% %%
% corre_value= correlate_signals(y,yestimate);%%calculate the correlation between the 2 signals, the higher the value the
% %more correlated or similar the signals are.
% str=sprintf('similarity= %f%%', corre_value);
% %%plot the results.
% figure;
% plot(yestimate,'r');
% hold on;
% plot(y,'b');
% grid on;
% xlabel('samples number');
% ylabel('System Response');
% dim = [0.2 0.6 0.3 0.3];
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% legend('Identified Model','Orginal Model');

%% system comparison and verification
A=[1 a];
B=[b];
%C = [1 r1 r2];
arx1red=idpoly(denomi,numerator);%,[1 r1 r2],1); 
zi=iddata(y(500+nk:1024+nk),u(500:1024),2); 
compare(zi,arx1red);

