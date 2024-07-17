
clc
clear
tic
format long
theta = [-90:.1:90];
K = length(theta);
M = 10;
i = 2;
f_0 = zeros(1,K);
f_1 = zeros(1,K);
f_2 = zeros(1,K);
f_3 = zeros(1,K);

for m = -M:M
    f_0  = exp(-1j *pi * m * sind(theta)) + f_0;
    f_1  = exp(-1i * m * pi/M) * exp(-1j *pi * m * sind(theta)) + f_1;
    f_2  = exp(-1i * m^2 * pi/M^2) * exp(-1j * pi * m * sind(theta)) + f_2;
    f_3  = exp(-1i * m^3 * pi/M^3) * exp(-1j * pi * m * sind(theta)) + f_3;
    i = i+1;
end

plot(theta,10*log10(f_0)); grid on; hold on;
plot(theta,10*log10(f_1)); grid on; hold on;
plot(theta,10*log10(f_2)); grid on; hold on;
plot(theta,10*log10(f_3)); grid on; hold on;