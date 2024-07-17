function round2 = Second_derivative (N,alpha) % sineta_sintheta = alpha

% numerator_fraction = pi * cos(eta) * (0.5 *N* sin(N*pi*alpha) .* sin(pi*alpha/2) - cos(pi*alpha/2) .* sin(N*pi*alpha/2).^2);

A = N^2/2 * (cos(N*pi*alpha/2).^2 .*  sin(pi*alpha/2).^2 - sin(N*pi*alpha/2).^2 .* sin(pi*alpha/2).^2);
B = -N * cos(pi * alpha/2) .* sin(N*pi*alpha) .* sin(pi*alpha/2) ;
C =  sin(N*pi*alpha/2).^2 .* (1/2 *  sin(pi*alpha/2).^2 + 3/2 *  cos(pi*alpha/2).^2 ) ;
numerator_fraction = pi.^2 * (A + B + C);

denumerator_fraction = N^2 * sin(pi*alpha/2).^4;

% P1 = -sin(eta + pi*alpha*(N-0.5)) * (1 + pi*(N-0.5)*cos(eta));
% P2 =  sin(eta + pi*alpha*(N+0.5)) * (1 + pi*(N+0.5)*cos(eta));
% P3 =  sin(eta - pi*alpha*(N+0.5)) * (1 - pi*(N+0.5)*cos(eta));
% P4 = -sin(eta - pi*alpha*(N-0.5)) * (1 - pi*(N-0.5)*cos(eta));
% 
% derivative_numerator_fraction_part1 = (N*pi/8) * (P1+P2+P3+P4); 
% derivative_numerator_fraction_part2 = (pi/2) * (- (pi/2 *cos(eta)-1)*sin(pi/2 * alpha - eta).* sin(N*pi*alpha/2).^2 + 2*cos(pi*alpha/2 - eta).*sin(N*pi*alpha/2).*cos(N*pi*alpha/2).*(N*pi/2)*cos(eta) - (pi/2 *cos(eta)+1)*sin(pi/2 * alpha + eta).* sin(N*pi*alpha/2).^2 + 2*cos(pi*alpha/2 + eta).*sin(N*pi*alpha/2).*cos(N*pi*alpha/2)*(N*pi/2)*cos(eta));
% derivative_numerator_fraction = derivative_numerator_fraction_part1 - derivative_numerator_fraction_part2;
% 
% derivative_denumerator_fraction = (3*pi/2) * N^2 * sin(pi*alpha/2).^2 * cos(eta) .* cos(pi*alpha/2);

% round2 = (derivative_numerator_fraction .* denumerator_fraction - derivative_denumerator_fraction .* numerator_fraction)./(denumerator_fraction .^ 2);

%%
% alpha_0 = find(alpha==0);
% if denumerator_fraction (alpha_0) == 0
%     numerator_fraction (alpha_0) = numerator_fraction (alpha_0-1);
%     denumerator_fraction (alpha_0) = denumerator_fraction (alpha_0-1);
% end
%%

round2 = numerator_fraction ./ denumerator_fraction;  
round2 (isnan(round2)) = 1;

end