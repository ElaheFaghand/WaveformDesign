function [q] = q_(h, x_hat, alpha, A_p)
%Q_ Summary of this function goes here
%   Detailed explanation goes here
    size_ = size(A_p);
    q=0;
    for  p = 1:size_(3)
        q = q + (A_p(:, :, p) * h(:,p));
    end
    q = q + x_hat * alpha / 2;
    % q = squeeze(q);
end

