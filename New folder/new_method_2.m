clc
clear
tic
format long
theta = [-90:.1:90];
M = 11;
gamma_L = sind(-30);
gamma_R = sind(30);
L_to_lambda = (gamma_R-gamma_L) * (M - 1)/4; 
z_to_lambda = linspace (-L_to_lambda ,L_to_lambda ,1000*M);
I_z = sinc(z_to_lambda);
abs_I_z = abs (I_z);
plot(z_to_lambda , I_z,"LineWidth",2), grid on , hold on
plot(z_to_lambda , abs_I_z,"LineWidth",1), grid on , hold on

pks = findpeaks(abs_I_z, z_to_lambda);
level = number(pks(1));
yline(level)

antenna_pos = linspace (-L_to_lambda ,L_to_lambda ,M);
Y_antenna = sinc (antenna_pos) ;
plot(antenna_pos, Y_antenna,'o',"LineWidth",6)

delta = 1e-3;
I_candidated = find(abs(abs_I_z - level)<delta);
z_candidated = z_to_lambda(I_candidated);
I_candidated = abs_I_z (I_candidated);
plot(z_candidated, I_candidated,'*',"LineWidth",6)

for i = 1: M
    [a(i),b(i)] = min(abs(antenna_pos(i) - z_candidated));
    Y_antenna_moving(i) = sinc(z_candidated(b(i)));
    sn = sign(Y_antenna_moving(i));
    displacement(i) = a(i) * sn
    plot(z_candidated(b(i)),Y_antenna_moving(i),'^',"LineWidth",6)
end 

R_initial = Y_antenna' * Y_antenna ;
P_initial = Beam_Pattern(M,theta,R_initial/sum(sum(R_initial)));
figure
plot(theta, 10*log10(P_initial)); grid on; axis([-80,80, -60, 5]); 

Y_antenna_moving = Y_antenna_moving.* exp(1i*displacement);
R = Y_antenna_moving' * Y_antenna_moving ;
P = Beam_Pattern(M,theta,R/sum(sum(R)));
hold on;
plot(theta, 10*log10(P)); grid on; axis([-90,90, -60, 5]);

% matrix = 1 * ones(M);
% for i = 1:M
%     matrix(i, i) = 1;
% end
% R_new = R.*matrix;
% P_new = Beam_Pattern(M,theta,R_new/sum(sum(R_new)));
% figure
% plot(theta, 10*log10(P_new)); grid on; axis([-90,90, -60, 5]);
toc