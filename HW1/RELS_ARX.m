%% EECE 574 Self-tuning
%% Author: Holguer A. Becerra
%% assignment 1.
%% Professor: Guy Dummont.
%% RELS Configurable

clear
clc
load data2017.mat


na=10;% system order Ay(t)
nb=10;% system order By(t)
nc=1;% system order Ce(t)
nk=1; % delay

u=data(:,1); %System Input
y=data(:,2); %System Output


%lenght of the SISO
m=length(y);
%% init variables for calculations
theta=zeros(1,na+nb+nc)'; % create a zero vector for the coefficients a and b, 2 for a 2 for b
%% ARX  model init values
init_a=[0.1276 0.123 0.1336 0.08193 0.1259 0.0301 0.03264 -0.03727 -0.06075 0.01514];
init_b=[1.5 -0.4166 -0.249 -0.1309 -0.1256 -0.01326 -0.1484 -0.07776 -0.1142 -0.1521];
for ho=1;na
    theta(ho)=init_a(ho);
end
for ho=1;nb
    theta(ho+na)=init_b(ho);
end


alpha=1e4; %%init factor
residual=zeros(1,m);
P=alpha*eye(na+nb+nc); % Covariance Matrix
lambda=1; % Forgetting Factor
lambda_inv=1/lambda; % inverse Forgetting Factor
error=0;
error_n=1;
%% for errors and parameters graphs
coeffs_historial=zeros(na+nb+nc,m);
coeffs_historial_error=zeros(na+nb+nc,m);
%% algorithm
for k=(na+nb+nc):m,%%sweeping out y
    %taking the new seeds for the coefficients estimation or observations,
    %2 samples behind
    xt=[];
    for order=1:na,
        xt=[xt -y(k-order)];
    end
    for order=1:nb,
        xt=[xt u(k-order)];
    end
    for order=1:nc,
        xt=[xt residual(k-order)];
    end
    %math the residual
    residual(k)=y(k)-xt*theta;
    x=xt';%x from xT
    %calculating P cvariance matrix
    P=lambda_inv*(P-(P*x*xt*P)/(lambda+xt*P*x));
    %taking the covariance matrix calculate and update the weights (Thetha, which are A(t), and B(t) coeffiecients)
    theta_prev=theta;
    theta=theta-P*x*(xt*theta-y(k));
               %____
                 %|
                 %.-> K(t+1)=P*x=P(t+1)*x(t+1)

    error=(error+residual(k)*residual(k));
    error=error/2;
    %error=sqrt(error/m);
    error_v(error_n)=error;
    error_n=error_n+1;

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
   % texto= [texto;sprintf('a%d',order)];
    hold on;
end
legend(texto);
xlabel('k')
ylabel('parameter')
texto=[];
subplot(2,2,2)
for order=1:(nb),
    plot(coeffs_historial(order+na,:));
   % texto= [texto;sprintf('b%d',order)];
    hold on;
end
legend(texto);
xlabel('k')
ylabel('parameter')
texto=[];
subplot(2,2,[3,4])
for order=1:(nc),
  plot(coeffs_historial(order+nb+na,:));
%  texto= [texto;sprintf('c%d',order)];
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
  %  texto= [texto;sprintf('a%d',order)];
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
  %  texto= [texto;sprintf('b%d',order)];
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
 %texto= [texto;sprintf('c%d',order)];
 pointText=['min c' num2str(order) '=' num2str(coeffs_historial_error(order+na+nb,m))];
    %text(m,coeffs_historial_error(order+na+nb,m)+0.1*order,pointText,'HorizontalAlignment','right');
 hold on;
end
legend(texto);
xlabel('k')
ylabel('Error')

%% 
figure;
plot(error_v);
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
%corre_value= correlate_signals(y,yestimate);%%calculate the correlation between the 2 signals, the higher the value the
%more correlated or similar the signals are.
%str=sprintf('similarity= %f%%', corre_value);
%%plot the results.
%figure;
%plot(yestimate,'r');
%hold on;
%plot(y,'b');
%grid on;
%xlabel('samples number');
%ylabel('System Response');
%dim = [0.2 0.6 0.3 0.3];
%annotation('textbox',dim,'String',str,'FitBoxToText','on');
%legend('Identified Model','Orginal Model');

%% system comparison and verification
A=[1 a];
B=[b];
C = [1 c];
arx1red=idpoly(A,B,C,1); 
zi=iddata(y(100+20+nk:501+nk),u(20+100:501),2); 
zi2=iddata(y(500+20+nk:1024+nk),u(20+500:1024),2); 
compare(zi,arx1red);
compare(zi2,arx1red);

