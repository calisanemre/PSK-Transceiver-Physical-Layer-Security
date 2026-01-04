classdef QPSKReceiver < handle
    properties
        os = 8;
        rolloff = 0.35;
        span = 10;
        L = 64;
        
        DATA_LEN = 500;
        DATA_SYMBOL_LEN;
        DATA_SEED = 123;
        tx_data_bits

        PLOT_LEN = 32;
        PLOT_SEED = 456;
        plots;

        DUMMY_BIT_LEN = 64;
        DUMMY_BIT_SEED = 789;
        
        ZEROS_LEN = 80;
        
        NOISE_SEED = 135;
        with_noise;
        alphaAN = 0.1;
        noise;
        noise_mode ='xor'; %  'phase_shift_an' 'gaussian_an'

        zc;
        zc_root = 25;
        zc_len  = 63;
        
        Nframe = 10;
        
        M = 4;
        phaseOffsetPSK = pi/4;

        rrc;

        fs;
        Ts;
        
        show_plots = false;
        before_time;
        after_time;
        after_cfo;
        after_phaseshift;
        after_nc;
    end
    
    methods
        function obj = QPSKReceiver()
            obj.fs = obj.os;
            obj.Ts = 1/obj.fs;
            obj.rrc = rcosdesign(obj.rolloff, obj.span, obj.os, 'sqrt');
            zcs = transpose(zadoffChuSeq(obj.zc_root, obj.zc_len));
            obj.zc = zcs(:);
            plot1_bits = obj.generateSeededBits(obj.PLOT_LEN, obj.M, obj.PLOT_SEED);
            obj.plots = obj.bits2psk(plot1_bits,obj.M, obj.phaseOffsetPSK);
            
        end
        
        function rxBits = receive(obj, rx)
            %% 1) Matched filter
            rx_filt = filter(obj.rrc,1,rx);
            rx_filt = rx_filt(obj.span*obj.os+1:end);
            rx_filt = rx_filt(:).';
            
           %% ===============================
            % MEYR–OEDER TIMING RECOVERY (ÖNCE)
            %% ===============================
            obj.before_time = rx_filt;
            [rx_sync, ~] = THAL_meyr_oeder_symbol_sync(obj.os, obj.L, rx_filt);
            obj.after_time = rx_sync;
            
            rx_sync = rx_sync(:);   % column vector
            
            %% ===============================
            % CFO CORRECTION (SONRA)
            %% ===============================
            
            % M-th power CFO estimator
            rx_cfo_est = rx_sync.^obj.M;   % QPSK → 4. kuvvet
            
            Nfft = 2^nextpow2(length(rx_cfo_est));
            RX = fftshift(fft(rx_cfo_est, Nfft));
            psd = abs(RX).^2;
            
            f = (-Nfft/2:Nfft/2-1)*(obj.fs/Nfft);
            
            [~, idx] = max(psd);
            f_hat = f(idx) / obj.M;
            
            n = (0:length(rx_sync)-1).';   % column vector
            rx_sync = rx_sync .* exp(-1j*2*pi*f_hat*n*obj.Ts);
            obj.after_cfo = rx_sync;
            
            %% ===============================
            % ZADFF-CHU
            %% ===============================

            corr_zc = abs(conv(rx_sync, conj(flipud(obj.zc))));

            % Tüm peak'leri bul (threshold olmadan)
            [all_pks, all_locs] = findpeaks( ...
                corr_zc, ...
                'MinPeakDistance', obj.zc_len, ...
                'SortStr', 'descend' ...  % En büyükten küçüğe sırala
            );
            
            if isempty(all_pks)
                error('No peaks found in ZC correlation');
            end
            
            % En büyük pikin %50'sinden büyük olanları al
            max_peak = all_pks(1);
            threshold = 0.5 * max_peak;  % %50 threshold (0.3-0.7 arası ayarlanabilir)
            
            valid_idx = all_pks >= threshold;
            pks = all_pks(valid_idx);
            locs = all_locs(valid_idx);
            
            % Lokasyona göre tekrar sırala (zaman sırasına göre)
            [locs, sort_idx] = sort(locs);
            pks = pks(sort_idx);

            frame_length = obj.zc_len + obj.PLOT_LEN + obj.DATA_SYMBOL_LEN;
            zc_starts = locs - obj.zc_len + 1;
            data_ends = zc_starts + frame_length - 1;
            
            % Filter: valid start AND valid end
            valid_mask = (zc_starts >= 1) & (data_ends <= length(rx_sync));
            
            pks = pks(valid_mask);
            locs = locs(valid_mask);

            if isempty(locs)
                error('No peaks with sufficient symbols (need %d symbols per frame)', frame_length);
            end

            if obj.show_plots
                figure;
                plot(corr_zc,'LineWidth',1.2);
                grid on;
                xlabel('Symbol index');
                ylabel('|Correlation|');
                title('Zadoff–Chu Correlation');
                hold on;
                xline(locs,'r--','Peak');
                legend('Correlation magnitude','Detected ZC peak');
            end
            
            n_peaks = length(pks);
            [zc_start, data_start] = deal(zeros(n_peaks,1));
            rx_data = zeros(n_peaks, obj.DATA_SYMBOL_LEN);
            rx_bits = zeros(n_peaks, log2(obj.M), obj.DATA_SYMBOL_LEN);            

            
            
            for i = 1:n_peaks
                zc_start(i) = locs(i) - obj.zc_len + 1;

                plot_start = zc_start(i) + obj.zc_len;
                plot_data = rx_sync(plot_start:plot_start + obj.PLOT_LEN -1);
                
                data_start(i) = zc_start(i) + obj.zc_len + obj.PLOT_LEN;
                rx_data_raw = rx_sync(data_start(i) : data_start(i) + obj.DATA_SYMBOL_LEN - 1).';
                
                a_hat = sum(plot_data .* conj(obj.plots)) / sum(abs(obj.plots).^2);
                
                data_start(i) = zc_start(i) ...
                   + obj.zc_len ...
                   + obj.PLOT_LEN; 
                
                rx_syncc = rx_sync/ a_hat;
                obj.after_phaseshift = rx_syncc;
                rx_data_corrected = rx_data_raw / a_hat;

                if obj.with_noise                  
                    switch lower(obj.noise_mode)
                        case 'phase_shift_an'
                            corrected = obj.substractArtificialNoise(rx_data_corrected, obj.alphaAN, obj.NOISE_SEED);
                            rx_data_corrected = corrected;
                            obj.after_nc = rx_data_corrected;
                            rx_data(i,:) = rx_data_corrected / sqrt(mean(abs(rx_data_corrected).^2));
                            rx_bits(i,:,:) = pskdemod(rx_data(i,:),obj.M,obj.phaseOffsetPSK,'OutputType','bit');
                        case 'gaussian_an'
                            corrected = obj.substractGaussianArtificialNoise(rx_data_corrected, obj.alphaAN, obj.NOISE_SEED);
                            rx_data_corrected = corrected;
                            obj.after_nc = rx_data_corrected;
                            rx_data(i,:) = rx_data_corrected / sqrt(mean(abs(rx_data_corrected).^2));
                            rx_bits(i,:,:) = pskdemod(rx_data(i,:),obj.M,obj.phaseOffsetPSK,'OutputType','bit');
                        case 'xor'
                            an = obj.generateSeededBits(obj.DATA_LEN/2, obj.M, obj.NOISE_SEED);
                            an2 = reshape(an, 2, []);
                            obj.after_nc = rx_data_corrected;
                            rx_data(i,:) = rx_data_corrected / sqrt(mean(abs(rx_data_corrected).^2));
                            noised = pskdemod(rx_data(i,:),obj.M,obj.phaseOffsetPSK,'OutputType','bit');
                            rx_bits(i,:,:) = xor(noised, an2);
                        otherwise
                        error('Unknown data_mode');
                    end
                else
                    obj.after_nc = rx_data_corrected;
                    rx_data(i,:) = rx_data_corrected / sqrt(mean(abs(rx_data_corrected).^2));
                    rx_bits(i,:,:) = pskdemod(rx_data(i,:),obj.M,obj.phaseOffsetPSK,'OutputType','bit');
                        
                end 
            end

            rxBits = permute(rx_bits(:,:,:), [1 3 2]);    
           
            
            if obj.show_plots
                fprintf('TX DATA symbols: %d | RX usable DATA symbols: %d\n', ...
                    obj.DATA_LEN, length(rx_data));
                figure;
                scatter(real(rx_data), imag(rx_data), 25,'filled');
                grid on; axis equal;
                xlabel('In-phase'); ylabel('Quadrature');
                title('RX DATA Constellation (After Phase Correction)');
            end
        
            
        end
        
        
        function [y, phi] = substractArtificialNoise(obj, x, alphaAN, seed)
            % x        : complex data symbols (row or column)
            % alphaAN  : phase noise variance (ör: 0.001, 0.01, 0.1)
            % seed     : RNG seed (optional)
            %
            % y        : phase-noised symbols
            % phi      : applied phase noise (rad)  <-- alıcı için
            
            if nargin > 2
                rng(seed);
            end
        
            % === 1) Generate phase noise (zero-mean) ===
            % alphaAN = Var{phi}
            phi = sqrt(alphaAN) * randn(size(x));   % radians
        
            % === 2) Apply phase-only rotation ===
            y = x .* exp(-1j * phi);
        end

        function y = substractGaussianArtificialNoise(obj, x, var, seed)

            rng(seed);
        
            an = sqrt(var/2) * (randn(size(x)) + 1j*randn(size(x)));
            y = x - an;
        
        end
    end
    methods (Static)

        function modSignal = bits2psk(bits, M, phaseOffset)
            % bits: 0-1 dizisi (column vector)
            % M: PSK modülasyon derecesi (2, 4, 8, ...)
            % phaseOffset: PSK phase offset, örn: pi/4
        
            k = log2(M); % Her semboldeki bit sayısı
            if mod(length(bits), k) ~= 0
                % Bit sayısı M ile bölünemiyorsa sıfır ekle
                bits = [bits; zeros(k - mod(length(bits), k), 1)];
            end
            
            % Bitleri gruplara ayır
            bitGroups = reshape(bits, k, []).';  % Her satır 1 sembol
            
            % Bin -> decimal
            symbols = bi2de(bitGroups, 'left-msb');  % Decimal sembol dizisi
            
            % PSK modülasyonu
            modSignal = pskmod(symbols, M, phaseOffset, 'gray');
        end


        function bits = generateSeededBits(dataLen, M, seed)
            % dataLen: üretilmek istenen sembol sayısı
            % M: PSK modülasyon derecesi
            % seed: random seed
        
            k = log2(M);  % Her semboldeki bit sayısı
            rng(seed);    % Seedli random üretim
            
            bits = randi([0 1], dataLen * k, 1);  % Random bit dizisi
            
        end

        function str = bits_linear2string(rxBits, originalBitCount)
            % BITS_LINEAR2STRING Decode Hamming(7,4) coded bits back to string
            %
            % INPUT:
            %   rxBits - received coded bit sequence (column vector)
            %   originalBitCount - original number of bits before padding (optional)
            % OUTPUT:
            %   str - decoded string
            
            % Hamming(7,4) parity check matrix
            H = [1 1 0 1 1 0 0;
                 1 0 1 1 0 1 0;
                 0 1 1 1 0 0 1];
            
            rxBits = reshape(rxBits, 7, []);  % 7-bit gruplar
            decodedBits = zeros(4, size(rxBits, 2));
            
            % Decode each block
            for i = 1:size(rxBits, 2)
                r = rxBits(:, i);
                s = mod(H * r, 2);  % syndrome
                syndrome_dec = bi2de(s', 'left-msb');
                
                if syndrome_dec ~= 0
                    % hatalı biti düzelt
                    r(syndrome_dec) = mod(r(syndrome_dec) + 1, 2);
                end
                decodedBits(:, i) = r(1:4); % sadece info bitler
            end
            
            % Bitleri tek sütun haline getir
            decodedBits = decodedBits(:);
            
            % ⭐ Eğer orijinal bit sayısı verilmişse, padding'i çıkar
            if nargin > 1 && ~isempty(originalBitCount)
                decodedBits = decodedBits(1:originalBitCount);
            end
            
            % 8-bit gruplara dönüştür ve karaktere çevir
            numCompleteBytes = floor(length(decodedBits) / 8);
            
            if numCompleteBytes > 0
                validBits = decodedBits(1:numCompleteBytes * 8);
                decodedBits = reshape(validBits, 8, []).';
                ascii_vals = bi2de(decodedBits, 'left-msb');
                str = char(ascii_vals).';
            else
                str = '';
            end
        end

               
    end
end