clf
clc
clear
tic

%% single lobe desired beampattern
start = -30;
stop = 30;
N = 10; % number of antennas

%%
theta_0 = (start + stop)/2;
theta = (-90:0.1:90)'*(pi/180); % radian
tt = (-90:0.1:90)'; % degree
L = length(theta);
for l = 1:L
    aT = exp(  1j * pi*(0:N-1)*sin(theta(l)) ).';
end

%%
a0 = exp(  1j * pi*(0:N-1)*sin(theta_0 *pi/180) ).';
beta0 = 1/(norm (a0)^2+norm (a0)^4);

IN = eye(N);
R_pre = beta0*(a0*a0'+IN);

%% Desired BeamPattern
pl1 = sin (start * pi/180 );
ph1 = sin (stop * pi/180);
P_d = (1.*(sin(theta)>pl1).*(sin(theta)<ph1));
theta_main = [start + 5 : 0.1 : stop - 5]' * (pi/180);

for i=1:length(P_d)
    if P_d(i)==0;
       P_d(i)=1e-1;
    end
end
plot(tt,10*log10(P_d)), grid on, hold on

%% Part_1
k=1;
D_k_while =1;
epsilon = 0.2;
eta = 0.1;
while max(D_k_while) > epsilon && k < 5
    BP_R_pre = Beam_Pattern (N,theta,R_pre);
    D_k = abs(P_d - BP_R_pre');
    [~,idx] = max (D_k);
    tt(idx)
    Ld_tetak = P_d(idx);
    aT_k = exp(  1j * pi*(0:N-1)*sin(theta(idx))).';
    beta_o = (Ld_tetak - (aT_k'*R_pre*aT_k))/norm (aT_k)^4;
    if beta_o >= -1/(aT_k'*R_pre^-1*aT_k);
        beta(k) = beta_o;
    else
         beta(k)= -1/(aT_k'*R_pre^-1*aT_k);
    end
    R(:,:,k) = R_pre + beta(k) * aT_k *aT_k';
    R_pre = R(:,:,k);
    P_sqp = Beam_Pattern (N,theta,R(:,:,k) );
   plot(tt,10*log10(abs(P_sqp)/max(P_sqp))), grid on, hold on

    BP_R_pre_while = Beam_Pattern (N,theta_main,R_pre);
    D_k_while = abs(ones(length(BP_R_pre_while),1) - BP_R_pre_while');
    k
    k=k+1;
    max(D_k_while)
end
figure
plot(tt,10*log10(abs(P_sqp)/max(P_sqp))), grid on, hold on

 %% Part 2
theta_side = [[-75 : 0.1 : start - 5],[ stop + 5 : 0.1 : 75]]' * (pi/180);
D_k_2 = 1;
k = 1;
while max(D_k_2) > eta && k < 10
        D_k_2 = abs(Beam_Pattern (N,theta_side,R_pre));
        [~,idx] = max (D_k_2);
        theta_side(idx) * 180/pi

        idx = find(abs(theta - theta_side(idx))<1e-5); % index in theta

        Ld_tetak = P_d(idx);
        aT_k = exp( 1j * pi*(0:N-1) * sin(theta(idx) )).';
        beta_o = (Ld_tetak - (aT_k'*R_pre*aT_k))/norm (aT_k)^4;
        if beta_o >= -1/(aT_k'*R_pre^-1*aT_k);
            beta(k) = beta_o;
        else
             beta(k)= -1/(aT_k' * R_pre^-1 * aT_k);
        end
        R(:,:,k) = R_pre + beta(k) * aT_k *aT_k';
        R_pre = R(:,:,k);
        P_sqp = Beam_Pattern (N,theta,R_pre );
        
        plot(tt,10*log10(abs(P_sqp)/max(P_sqp))), grid on, hold on
        k = k+1;
end
%%
function Beam_Pattern = Beam_Pattern (N,gamma,R)
Beam_Pattern = zeros(1, length(gamma));
    for l = 1:length(gamma)
        et = exp(1j*pi*(0:N-1)*(gamma(l))).';
        Beam_Pattern(l) = et'*R*et;
    end
end
