

clear               %??
 
 a(1)=1;b(1)=0;c(1)=1;d(1)=0;u(1)=d(1);z(1)=0;z(2)=0; %???
% for i=2:1200          %??m??u(i)         
%     a(i)=xor(c(i-1),d(i-1));
%     b(i)=a(i-1);
%     c(i)=b(i-1);
%     d(i)=c(i-1);
%     u(i)=d(i);
% end
% u; %???‘?’???????????m??
%  
% v=randn(1200,1); %?????????
load data2017.mat

u=data(:,1); %System Input
z=data(:,2); %System Output
m=length(u);
v=zeros(1,m);
V=0;   %??????
for i=1:m
    V=V+v(i)*v(i);
end
V1=V/m;   

 
o1=0.001*ones(6,1);p0=eye(6,6);         %???
zf(1)=0.1;zf(2)=0.1;vf(2)=0.1;vf(1)=0.1;uf(2)=0.1;uf(1)=0.1;
   %???????????
for k=3:m
      h=[-z(k-1);-z(k-2);u(k-1);u(k-2);v(k-1);v(k-2)];
      hf=h;
     K=p0*hf*inv(hf'*p0*hf+1); 
     p=[eye(6,6)-K*hf']*p0; 
     v(k)=z(k)-h'*o1;
     o=o1+K*v(k) ;
       p0=p;
    o1=o;
    a1(k)=o(1);
    a2(k)=o(2);
    b1(k)=o(3);
    b2(k)=o(4);
    d1(k)=o(5);
    d2(k)=o(6);
    e1(k)=abs(a1(k)+1.2);
    e2(k)=abs(a2(k)-0.6);
    e3(k)=abs(b1(k)-1.0);
    e4(k)=abs(b2(k)-0.5);
    e5(k)=abs(d1(k)+1.0);
    e6(k)=abs(d2(k)-0.2);
     zf(k)=z(k)-d1(k)*zf(k-1)-d2(k)*zf(k-2);
 
     uf(k)=u(k)-d1(k)*uf(k-1)-d2(k)*uf(k-2);   
     
     vf(k)=v(k)-d1(k)*vf(k-1)-d2(k)*vf(k-2);  
   hf=[-zf(k-1);-zf(k-2);uf(k-1);uf(k-2);vf(k-1);vf(k-2)];    
 end 
 o1  %???‘?’?????????????
 V1
  %??
 subplot(4,1,1)
k=1:m;
plot(k,a1,'k:',k,a2,'b',k,b1,'r',k,b2,'m:',k,d1,'g',k,d2,'k');
xlabel('k')
ylabel('parameter')
legend('a1=-1.2,','a2=0.6','b1=1.0','b2=0.5','d1=-1.0','d2=0.2'); %???
title('The parameter idendification of the RML');

subplot(4,1,2)
k=1:m;
plot(k,e1,'k',k,e2,'b',k,e3,'r',k,e4,'m',k,e5,'g',k,e6,'k');
xlabel('k')
ylabel('error')
%title('????')

subplot(4,1,3)
k=1:m;
plot(k,u);
xlabel('k')
ylabel('input')
%title('??????')

subplot(4,1,4)
k=1:m;
plot(k,v);
xlabel('k')
ylabel('random noise')
%title('?????????')

A=[1 o(1) o(2)];
B=[o(3) o(4)];
C = [1 o(5) o(6)];
arx1red=idpoly(A,B,C,1); 
zi=iddata(z(100+21:502),u(20+100:501),1); 
%zi2=iddata(y(500+21:1025),u(20+500:1024),2); 
compare(zi,arx1red);