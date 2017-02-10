%% EECE 574 Self-tuning
%% Author: Holguer A. Becerra
%% assignment 1.
%% Professor: Guy Dummont.
%% RML Configurable

clear
clc
load data2017.mat


na=3;% system order Ay(t)
nb=2;% system order By(t)
nc=1;% system order Ce(t)
nk=1; % delay

u=data(:,1); %System Input
y=data(:,2); %System Output

%lenght of the SISO
m=length(u);

%% init variables for calculations
theta=zeros(1,na+nb+nc)'; % create a zero vector for the coefficients a and b, 2 for a 2 for b
%% transfer function estimations to optimize
init_a=[-0.5784 0.002871 -0.01826];
init_b=[1.5 -1.475];
theta(1)=init_a(1);
theta(2)=init_a(2);
theta(3)=init_a(3);
theta(4)=init_b(1);
theta(5)=init_b(2);

alpha=1e4; %%init factor
residual=zeros(1,m);
P=alpha*eye(na+nb+nc); % Covariance Matrix
lambda=1; % Forgetting Factor
lambda_inv=1/lambda; % inverse Forgetting Factor
n_residuals_last=0.1;
n_residuals=0;
yhat=zeros(1,m);
uhat=zeros(1,m);
rhat=zeros(1,m);
d=[];
for init=1:m,
    yhat(init)=0.1;
    uhat(init)=0.1;
    rhat(init)=0.1;
    residual(init)=0.1;
end
%% for errors and parameters graphs
coeffs_historial=zeros(na+nb+nc,m);
coeffs_historial_error=zeros(na+nb+nc,m);
%% algorithm
for k=(na+nb+nc):m,%%sweeping out y
    %taking the new seeds for the coefficients estimation or observations,
    %2 samples behind
    %% collect u and y to form x
    xt=[];
    %n_residuals=[];
    for order=1:na,
        xt=[xt -y(k-order)];
    end
    for order=1:nb,
        xt=[xt u(k-order)];
    end

    for order=1:nc,
        xt=[xt residual(k-order)];
        %store the previous residuals in d
        d=[d residual(k-order)];
    end

    x=xt';%this would be x from x Transpose
    % calculate the residula 
    residual(k)=y(k)-xt*theta;

    %% compute yhat and uhat and residual hat
    % using d progress
    for order=1:nc,
     yhat(k)=y(k)-d(order)*yhat(k-order);
    end
    for order=1:nc,
     uhat(k)=u(k)-d(order)*uhat(k-order);
    end
    for order=1:nc,
     rhat(k)=residual(k)-d(order)*rhat(k-order);
    end
    xtf=[];
    %% generate xf 
    for order=1:na,
        xtf=[xtf -yhat(k-order)];
    end
    for order=1:nb,
        xtf=[xtf uhat(k-order)];
    end
    for order=1:nc,
        xtf=[xtf rhat(k-order)];
    end

    %% RML calculation
    xf=xtf';
    P=P-((P*xf*xtf*P)*inv(1+xtf*P*xf));
    K=P*xf*inv(1+xtf*P*xf);
    theta_prev=theta;
    theta=theta+K*(residual(k));

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

     for order=1:nc,
         coeffs_historial(order+nb+na,k)=theta(order+nb+na);
         coeffs_historial_error(order+nb+na,k)=abs(coeffs_historial(order+nb+na,k)-theta_prev(order+nb+na));
     end
end

%% Graph Parameters
texto=[];
figure;
subplot(2,2,1)
for order=1:(na),
    plot(coeffs_historial(order,:));
    texto= [texto;sprintf('a%d',order)];
    hold on;
end
legend(texto);
xlabel('k')
ylabel('parameter')
texto=[];
subplot(2,2,2)
for order=1:(nb),
    plot(coeffs_historial(order+na,:));
    texto= [texto;sprintf('b%d',order)];
    hold on;
end
legend(texto);
xlabel('k')
ylabel('parameter')
texto=[];
subplot(2,2,[3,4])
for order=1:(nc),
  plot(coeffs_historial(order+nb+na,:));
  texto= [texto;sprintf('c%d',order)];
  hold on;
end

legend(texto);
xlabel('k')
ylabel('parameter')


%% GRAPH Errors
texto=[];
figure;
subplot(2,2,1)
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
subplot(2,2,2)
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
texto=[];
subplot(2,2,[3,4])
for order=1:(nc),
 plot(coeffs_historial_error(order+nb+na,:));
 texto= [texto;sprintf('c%d',order)];
 pointText=['min c' num2str(order) '=' num2str(coeffs_historial_error(order+na+nb,m))];
    %text(m,coeffs_historial_error(order+na+nb,m)+0.1*order,pointText,'HorizontalAlignment','right');
 hold on;
end
legend(texto);
xlabel('k')
ylabel('Error')

%% 
%figure;
%plot(error_v);
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
C = [1 c];
arx1red=idpoly(A,B,C,1); 
zi=iddata(y(100+20+nk:501+nk),u(20+100:501),2); 
zi2=iddata(y(500+20+nk:1024+nk),u(20+500:1024),2); 
compare(zi,arx1red);
compare(zi2,arx1red);
