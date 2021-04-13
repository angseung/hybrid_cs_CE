clc; clear; close all;

K = 128; % number of OFDM symbols
M = 4;   % swipping level
N = 20;  % number of bases in a basis set
Nr = 8;  % number of UE antennas in an RF chain

seq = gen_zadoffchu(K + 1, M);

%% Check seq orthogonality...
% a = seq(1, 1 : fft_size);
% b = seq(2, 1 : fft_size);
% c = seq(3, 1 : fft_size);
% d = seq(4, 1 : fft_size);
% 
% abs(a * b')
% abs(a * c')
% abs(a * d')
% 
% abs(b * c')
% abs(b * d')
% abs(c * d')

seq = seq(:, 1 : K);

%% Calculate W

W_idt = [eye(M), zeros(M, Nr - M)];
U = eye(M);

W = U * W_idt;
