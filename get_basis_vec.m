function [a_theta] = get_basis_vec(Nr, theta)
    % return : (Nr, 1) shape matrix
    exp_theta = 0 : Nr - 1;
    a_theta = exp(1j * pi * cos(theta) .* exp_theta).';
    
end