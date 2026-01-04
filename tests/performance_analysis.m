%% QPSK System Performance Analysis
% This script tests the performance of a QPSK system built on the ITU THAL 
% framework under the influence of Artificial Noise (AN).

clc; clear; close all;

%% 1. Parameter Configuration
EbNo = 30;                  % Channel noise level (dB)
alphaAN_vec = 0.1:0.3:1.0;  % Artificial Noise variance range (Sweep)

% Initialize Objects
tx = QPSKTransmitter();
rx_target = QPSKReceiver(); % Authorized receiver (Possesses key to remove AN)
rx_hidden = QPSKReceiver(); % Unauthorized/Hidden receiver (No key to remove AN)

% Noise Mode and Activation
tx.with_noise = true;
tx.noise_mode = 'gaussian_an'; % Options: 'gaussian_an', 'phase_shift_an', 'xor'

rx_target.noise_mode  = tx.noise_mode;
rx_hidden.noise_mode  = tx.noise_mode;

rx_target.with_noise = true;   % Noise cancellation enabled
rx_hidden.with_noise = false;  % Noise cancellation disabled (Eavesdropper scenario)

%% 2. Simulation Loop
BER_target = zeros(size(alphaAN_vec));
BER_hidden = zeros(size(alphaAN_vec));

fprintf('Starting Simulation (Eb/No = %d dB)...\n', EbNo);

for k = 1:length(alphaAN_vec)
    % Update AlphaAN
    current_alpha = alphaAN_vec(k);
    tx.alphaAN = current_alpha;
    rx_target.alphaAN = current_alpha;
    rx_hidden.alphaAN = current_alpha;

    %% --- TRANSMIT ---
    tx_signal = tx.transmit();

    %% --- CHANNEL EFFECTS ---
    % 1. Carrier Frequency Offset (CFO) Injection
    fs = tx.os * 1e6;
    Delta_f = 0.001 * fs; % Example frequency shift
    Ts = 1/fs;
    n = (0:length(tx_signal)-1).';
    tx_cfo = tx_signal .* exp(1j*2*pi*Delta_f*n*Ts);

    % 2. Add AWGN (Channel Noise)
    rx_channel = awgn(tx_cfo, EbNo, 'measured');

    %% --- RECEIVE ---
    % Synchronize parameters
    rx_target.DATA_LEN = tx.DATA_LEN;
    rx_hidden.DATA_LEN = tx.DATA_LEN;
    rx_target.DATA_SYMBOL_LEN = tx.DATA_SYMBOL_LEN;
    rx_hidden.DATA_SYMBOL_LEN = tx.DATA_SYMBOL_LEN;

    % Execute Receivers
    rx_target_bits = rx_target.receive(rx_channel);
    rx_hidden_bits = rx_hidden.receive(rx_channel);
    
    [N_frames, ~, ~] = size(rx_target_bits);

    %% --- BIT ERROR RATE (BER) CALCULATION ---
    % Reshape matrices into vectors for comparison
    rx_t_vec = reshape(permute(rx_target_bits,[3 2 1]), [], 1);
    rx_h_vec = reshape(permute(rx_hidden_bits,[3 2 1]), [], 1);
    
    % Replicate transmitted bits based on the number of frames received
    tx_bits_rep = repmat(tx.tx_data_bits(:), N_frames, 1);

    BER_target(k) = mean(rx_t_vec ~= tx_bits_rep);
    BER_hidden(k) = mean(rx_h_vec ~= tx_bits_rep);
    
    fprintf('Alpha: %.2f | Target BER: %.4e | Hidden BER: %.4e\n', ...
            current_alpha, BER_target(k), BER_hidden(k));
end

%% 3. Visualization
figure('Name', 'QPSK Artificial Noise Analysis');

% Linear Plot
subplot(2,1,1);
plot(alphaAN_vec, BER_target, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
plot(alphaAN_vec, BER_hidden, 'r-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
grid on;
xlabel('\alpha_{AN} (Noise Power)'); ylabel('BER');
legend('Target Receiver (Keyed)', 'Hidden Receiver (No Key)');
title('BER vs. Artificial Noise Power (Linear Scale)');

% Log-Log Plot
subplot(2,1,2);
loglog(alphaAN_vec, max(BER_target, 1e-6), 'b-o', 'LineWidth', 1.5); hold on;
loglog(alphaAN_vec, max(BER_hidden, 1e-6), 'r-s', 'LineWidth', 1.5);
grid on;
xlabel('\alpha_{AN}'); ylabel('BER (Log Scale)');
title('BER vs. Artificial Noise Power (Logarithmic Scale)');