function R_shift = shift_R_construction (of,N,R)
    q = zeros(N);
        for i=1:N
            q(i,i)= exp(1i*(pi)*(i-1)*sin(of));
        end
    R_shift=q*R*q';
end