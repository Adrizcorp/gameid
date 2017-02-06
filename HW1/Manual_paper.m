N=length(data(:,2));
h=1;%time sampling
N1=N/2;
y=data(:,2);
u=data(:,1);
originales= iddata(y(1:N),u(1:N),h); 
zi=iddata(y(1:N1),u(1:N1),h); 
zv=iddata(y(N1+1:N),u(N1+1:N),h); 
x=[1:N1];
zi=detrend(zi); 
zv=detrend(zv); 
figure;
plot(zi);

yi=zi.OutputData; ui=zi.InputData; 

%spectrum(ui,yi); 
%S=spectrum(ui,yi); 
%delay estimation
nk=delayest(zi);
Wn=[0 0.4];
zif=idfilt(zi,6,Wn); 
figure;
plot(zif);

ir=cra(zi) 
sr=cumsum(ir) 
plot(sr)

na=1; nb=1; nc=1; nd=1; nf=1; nk=1; 
orders=[na nb nc nd nf nk]; 
pem1=pem(zi,orders) 
present(pem1); 

arx1=arx(zi,[10 9 1]); 
[A,B,C,D]=ssdata(arx1); 
[Ab,Bb,Cb,M,T]=dbalreal(A,B,C); 
[Ab,Bb,Cb,Db]=dmodred(Ab,Bb,Cb,D,4:10); 
[b,a]=ss2tf(Ab,Bb,Cb,Db); 
arx1red=idpoly(a,b); 
NN=[2 1 1;3 1 1;3 2 1;4 3 1;5 4 1]; 
V=arxstruc(zi,zv,NN); 
selstruc(V);
