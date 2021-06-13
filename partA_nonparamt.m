%-----------------Part A - Non - Paramteric Methods ---------------%
%Author: Awab Syed 
%Date 30 Jan 2021 
%Module: ACS 6129 - System Identification 

%-----------------------------PRBS Signal-----------------------------%
sim('prbs7.mdl')
% Calculate the autocorrelation of generated PRBS7 signal
y7 = double(y7);
y7(y7 == 0) = -1;           % replace 0 with -1
c = xcorr(y7)/length(y7);   % calculate correlation

% Plot autocorrelation
t = -126:1:126;
figure(1)
plot(t,c);
grid on
xlabel('Time Lag [s]');
ylabel('Autocorrelation');
title('Autocorrelation of PRBS Signal')

%---------------------Averaging process 50 Realizations------------%
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

figure(2)
hold on
plot(t,Ctemp)
plot(t,C,'g','LineWidth',2) 
hold off
xlabel('time lag');
ylabel('autocorrelation');
title ('Average Correlation of 50 Realisations')

%---------------------Step Testing---------------------------------%
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
figure(3)
plot(t,y,'k',t,avg_yn,'k:',t,sd_po,'r--',t,sd_ne,'r--');
legend('true','estimate')
xlabel('t');
ylabel('y(t)');

%-----------------------End of Part A-----------------------------%