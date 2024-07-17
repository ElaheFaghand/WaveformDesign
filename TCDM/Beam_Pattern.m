function Beam_Pattern = Beam_Pattern (N,gamma,R)
Beam_Pattern = zeros(1, length(gamma));
    for l = 1:length(gamma)
        et = exp(1j*pi*(0:N-1)*(gamma(l))).';
        Beam_Pattern(l) = et'*R*et;
    end
end