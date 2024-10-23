% Eq.36
function round2 = Second_derivative (N,alpha) % sineta_sintheta = alpha

A = N^2/2 * (cos(N*pi*alpha/2).^2 .*  sin(pi*alpha/2).^2 - sin(N*pi*alpha/2).^2 .* sin(pi*alpha/2).^2);
B = -N * cos(pi * alpha/2) .* sin(N*pi*alpha) .* sin(pi*alpha/2) ;
C =  sin(N*pi*alpha/2).^2 .* (1/2 *  sin(pi*alpha/2).^2 + 3/2 *  cos(pi*alpha/2).^2 ) ;
numerator_fraction = pi.^2 * (A + B + C);

denumerator_fraction = N^2 * sin(pi*alpha/2).^4;

round2 = numerator_fraction ./ denumerator_fraction;  
round2 (isnan(round2)) = 1;

end
