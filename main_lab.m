% Generate the PRBS3 signal
sim('prbs3.mdl') % Run 'prbs3.mdl'

% Plot the PRBS3 signal
t = 0:0.01:14.99;
sig = repelem(y3,100);
figure()
plot(t,sig)
grid on
xlabel('Time [s]')
ylabel('Amplitude')
title('PRBS3 Signal')
axis([0 15 -0.1 1.1])
%%
%--------------------Question 2-----------------------------%
% Calculate the autocorrelation of generated PRBS7 signal
y7 = double(y7);
y7(y7 == 0) = -1;           % replace 0 with -1
c = xcorr(y7)/length(y7);   % calculate correlation

% Plot autocorrelation
t = -126:1:126;
figure(2)
plot(t,c);
grid on
xlabel('Time Lag [s]');
ylabel('Autocorrelation');
title('Autocorrelation of PRBS Signal')
%axis([-126 126 -0.1 1.1])
%%
%-------------------Question 3-------------------------------%
% Calculates and plots the average autocorrelation
n = 7;
Ny = 2^n - 1;       % size of PRBS
Nc = 2*Ny + 1;      % size of the correlation array
ns = 50;            % number of realizations
Ctemp = zeros(Nc,ns);

for i=1:ns
    
    % Generate a random initial condition
    initial_state = rand(1,n)<=0.5;
    
    % Run prbs7.mdl
    sim('prbs7_avg.slx')
    y7 = double(y7);
    y7(y7 == 0) = -1;

    % Calculate autocorrelation
    c1 = xcorr(y7)/length(y7);

    % Calculate the average autocorrelation
    Ctemp(:,i) = c1;
end

C = mean(Ctemp,2);
t = -127:127;

figure()
hold on
plot(t,Ctemp)
plot(t,C,'g','LineWidth',2) 
hold off
xlabel('time lag');
ylabel('autocorrelation');
title ('Average Correlation of 50 Realisations')
%%
%-----------------------Question 4-----------------------------%
% Calculates and plots the power spectral density
sampling_frequency = 1000;
t = 0:1/sampling_frequency:10; % 10 seconds of data @ 1KHz
data1 = sin(t*50*2*pi); % 50 Hz
data2 = sin(t*100*2*pi); % 100 Hz
data = data1 + data2;

% power spectral density using Welch s method
[P,f] = pwelch(data,256,250,length(data),sampling_frequency);

figure(3)
subplot(211)
plot(t,data);
axis([0,.5,-2,2]);
xlabel('Time (s)');
ylabel('Amplitude');
title('sum of sinusoidal signals')

subplot(212)
plot(f,P);
title('power spectral density')
xlabel('Frequency (Hz)')
ylabel('Power')
%%
% Single run 
% System
num_sys = [0 1];
den_sys = [1 -.9];

% Observation noise
num_noise = 1;
den_noise = [1 -.9];

% Time properties
t0 = 0;
tfinal = 10;
sampling_time = 0.1;
t = t0:sampling_time:tfinal;

% Input and noise signals
n = length(t);
u = ones(n,1);
sigma = 1;
% e = sigma*randn(n,1);
% yn = y + filter(num_noise,den_noise,e);

% Simulation
y = filter(num_sys,den_sys,u);
for i = 1:100
    e = sigma*randn(n,1);
    yn(:,i) = y + filter(num_noise,den_noise,e);
end
%Average and Confidence level of 95%
avg_yn = mean(yn,2);
sd_po = avg_yn+1.96*std(yn,0,2)/sqrt(n);
sd_ne = avg_yn-1.96*std(yn,0,2)/sqrt(n);
% Plot
plot(t,y,'k',t,avg_yn,'k:',t,sd_po,'r--',t,sd_ne,'r--');
legend('true','estimate')
xlabel('t');
ylabel('y(t)');
%% Parametric Method 
n = 10;
% system
num_sys = [0 1];
den_sys = [1 -0.9];
% white binary input
u = sign(randn(n,1));

% output
y = filter(num_sys,den_sys,u);

e = 0;
yn = y+e;

% construct Y~ and Phi
Y = yn(2:end);
Phi = [-yn(1:end-1) u(1:end-1)];

% find beta
beta1 = inv(Phi'*Phi)*Phi'*Y
beta2 = Phi\Y
beta3 = pinv(Phi)*Y

%%
% ARX model 
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
plot(res_cor);
title('Autocorrelation with = 100')
xlabel('length')
ylabel('autocorrelation')
var_uncorrnoise = std(res)^2

%% Output Error Model 
n = 100;
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
plot(res_cor);
title('Autocorrelation with N = 100')
xlabel('length')
ylabel('autocorrelation')
var_uncorrnoise = std(res)^2
%%
%Histogram (use beta2)
m = 100;
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
figure(2)
subplot(2,1,1)
plot(Ahat);
xlabel('iteration');
ylabel('ahat');
subplot(2,1,2);
plot(Bhat);
xlabel('iteration');
ylabel('bhat')




