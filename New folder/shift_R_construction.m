function R_shift = shift_R_construction (of,N,R)
    A = zeros(N);
        for i=1:N
            A(i,i)= exp(1i*(pi)*(i-1)*(sin(of*pi/180)));
        end
    R_shift=A*R*A';
end