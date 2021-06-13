%-----------------Part B - Paramteric Methods ---------------%
%Author: Awab Syed 
%Date 30 Jan 2021 
%Module: ACS 6129 - System Identification 

%-------------------------------ARX Model--------------------------%
n = 100;

num_sys = [0 1];
den_sys = [1 -0.9];

u = sign(randn(n,1));

num_noise = 1;
den_noise = [1 -0.9];
sigma = 1;
e = sigma*randn(n,1);

y = filter(num_sys,den_sys,u);
yn = y + filter(num_noise,den_noise,e);

Y = yn(2:end);
Phi = [-yn(1:end-1) u(1:end-1)];

beta1 = inv(Phi'*Phi)*Phi'*Y
beta2 = Phi\Y;
beta3 = pinv(Phi)*Y;

Yhat = Phi*beta1;

res = Y-Yhat;
res_cor = xcorr(res)/length(res);
figure (1)
plot(res_cor);
title('Autocorrelation with = 100')
xlabel('length')
ylabel('autocorrelation')
var_uncorrnoise = std(res)^2
%-----------------------------Output Error Model----------------------%
n = 100;   %n = number of points/samples - change this to 10 if n=10
u = sign(randn(n,1));

num_sys = [0 1];
den_sys = [1 -.9];

sigma = 1;
e = sigma*randn(n,1);

y = filter(num_sys,den_sys,u);
yn = y+e;

Y = yn(2:end);
Phi = [-yn(1:end-1) u(1:end-1)];

beta1 = inv(Phi'*Phi)*Phi'*Y;
beta2 = Phi\Y;
beta3 = pinv(Phi)*Y;

Yhat = Phi*beta1;

res = Y-Yhat;
res_cor = xcorr(res)/length(res);
figure(2)
plot(res_cor);
title('Autocorrelation with N = 100') %Change N = 10 if n = 10 above 
xlabel('length')
ylabel('autocorrelation')
var_uncorrnoise = std(res)^2
%--------------------------------Histogram-----------------------%
m = 100;   %Number of runs 
Bhat = zeros(m,1);
Ahat = zeros(m,1);
beta = beta2

for i=1:m
    ahat = beta(1);
    Ahat(i) = beta(1);
    Bhat(i) = beta(2);
    
    den_Fhat = [1 ahat];
    num_Fhat = [1];
    
    yf = filter(num_Fhat,den_Fhat,yn);
    Yf = yf(2:end);
    
    uf = filter(num_Fhat,den_Fhat,u);
    
    Phif = [-yf(1:end-1),uf(1:end-1)];
    
    beta = inv(Phif'*Phif)*Phif'*Yf;
end
figure(3)
subplot(2,1,1)
plot(Ahat);
xlabel('iteration');
ylabel('ahat');
subplot(2,1,2);
plot(Bhat);
xlabel('iteration');
ylabel('bhat')

%---------------------------------End of Part B---------------------%