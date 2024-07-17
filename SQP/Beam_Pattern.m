function Beam_Pattern = Beam_Pattern (N,theta,R)
Beam_Pattern = zeros(1, length(theta));
    for l = 1:length(theta)
        et = exp(1j*pi*(0:N-1)*sin(theta(l))).';
        Beam_Pattern(l) = et'*R*et;
    end
end