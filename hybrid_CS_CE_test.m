clc; clear; close all;

K = 128; % number of OFDM symbols
M = 4;   % swipping level
N = 20;  % number of bases in a basis set
Nr = 8;  % number of UE antennas in an RF chain
path = 7; % number of paths
scatter = 10;
mod_type = 2;
snr = 10; % in dB scale

seq = gen_zadoffchu(K + 1, M);
seq = seq(:, 1 : K);
% seq : (subcarrier k, symbols) shape...
% then P = diag(seq(i, K))...

%% Calculate W
FIRST_OPT = 0;
SECOND_OPT = ~FIRST_OPT;

U = randU(M); % Any unitary matrix 'U' could be...

if (FIRST_OPT)
    W_idt = [eye(M), zeros(M, Nr - M)];
    
elseif (SECOND_OPT)
    b1 = 0;
    bm = b1;
    W_idt = zeros(M, Nr);
    
    for i = 1 : M
        W_idt(i, :) = get_basis_vec(Nr, bm)';
        bm = bm + (2 * pi * i) / M;
    end
    
end

W = U * W_idt;

%% Calculate true channel H

% Set default params
model = SCM();
model.n_path = path;
model.n_mray = scatter;
cp_len = 0;
data_len = K * mod_type;
tx_ant = 64;
rx_node = 1; % number of UEs

% Channel setting
N_d = rx_node;
N_s = rx_node;
A = [1];

% Set antenna shape
[~, ~] = model.ant(Nr, tx_ant);
N_tx = model.Ntx;
N_rx = model.Nrx;

% Get channel weight
r_H = zeros(path, K + cp_len, N_rx * N_d, N_tx);
% t_H = zeros(path, N_rx * N_d, N_tx);

for d = 1 : N_d
    temp = model.FD_channel(K + cp_len);
    r_H(:,:,1+(d-1)*N_rx:d*N_rx,:) = temp;
end

h = squeeze(r_H(:,1,:,:));
H = fft(h, K, 1);
H = squeeze(H(:, :, 1)).';

%% Calculate received signal Z
Z = zeros(M, K);

for i = 1 : M
    P = diag(seq(i, :));
    Z(i, :) = W(i, :) * awgn_noise(H * P, snr);
end

Y = Z / P;

%% Calculate X_bar for a basis set
X_bar = pi * (-1 + (2 .* (0 : N - 1)) / N );

%% Calculate estimated basis set x_hat
candidate_x = zeros(1, N);

for i = 1 : N
    curr_psi = W * get_basis_vec(Nr, X_bar(i));
    candidate_x(1, i) = ...
        (curr_psi' * Y * Y' * curr_psi) ./ (curr_psi' * curr_psi);
end

[x_hat_1, ~] = max(abs(candidate_x));

x_hat_m = x_hat_1 + 2 * pi .* (1 : M - 1) ./ M;
x_hat_m = [x_hat_1, x_hat_m];

%% Calculate sparse vector 'S'
basis = zeros(M, M);
for i = 1 : M
    basis(:, i) = W * get_basis_vec(Nr, x_hat_m(i));
end

S = basis' * Y;

%% Decise final estimated channel coefficient
A = zeros(Nr, M);
for i = 1 : M
    A(:, i) = get_basis_vec(Nr, x_hat_m(i));
end

H_est = A * S;
