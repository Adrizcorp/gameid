%************************************************************************%
clear;
clc;
load('data2017');
u = data(:,1); y = data(:,2);
u = detrend(u);
y = detrend(y);
data = iddata(y,u,1,'InputName','u(t)','OutputName','y(t)');
advice(data)
isnlarx(data,[1,2,3]);
ir=cra(data) 
sr=cumsum(ir) 
figure;
plot(sr);

%sys = arx(data, [2 2 2]);
dt = 1;
N = length(u);



%% create an ARX model
estimate = data(1:700);
verif = data(701:1025);
v = arxstruc(estimate,verif,struc(1:10, 1:10, 1:10));
% give me advice about the order I could take
order = selstruc(v,'aic');
% base on that start from 1 to the ARX suggestion order for ARMAX
na = 1:order(1);
nc = 1:order(3);
nk = 0:2;
L = length(na) * length(nc) * length(nk);
models = cell(1,L);
ct = 1;
for i = 1:length(na)
    na_ = na(i);
    nb_ = na_;
    for j = 1:length(nc)
        nc_ = nc(j);
        for k = 1:length(nk)
            nk_ = nk(k); 
            models{ct} = armax(data,[na_ nb_ nc_ nk_]);
            ct = ct+1;
        end
    end
end  

%% 
% Stack the estimated models and compare their simulated responses to estimation
% data |z|. 
models = stack(1,models{:});
compare(data,models)  