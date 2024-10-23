% Eq.35
function round1 = First_derivative (N,alpha) % sineta_sintheta = alpha
numerator_fraction = pi * sin (N * pi * alpha /2) .* (-N* cos(N*pi*alpha/2) .* ...
    sin(pi*alpha/2) + cos(pi*alpha/2) .* sin(N*pi*alpha/2));

denumerator_fraction = N^2 .* sin(pi*alpha/2).^3;

round1 = numerator_fraction./denumerator_fraction;

round1 (isnan(round1)) = 1;
end
