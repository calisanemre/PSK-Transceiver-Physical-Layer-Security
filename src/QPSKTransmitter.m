classdef QPSKTransmitter < handle
    properties
        os = 8;
        rolloff = 0.35;
        span = 10;
        
        DATA_LEN = 500;
        DATA_SYMBOL_LEN;
        DATA_SEED = 123;
        ORG_DATA_LEN;
        tx_data_bits;
        data_mode = 'analysis'; % 'visual'
        noise_mode = 'xor'; % 'phase_shift_an' 'gaussian_an'
        PLOT_LEN = 32; 
        PLOT_SEED = 456;
        plot;

        DUMMY_BIT_LEN = 64;
        DUMMY_BIT_SEED = 789;
        
        ZEROS_LEN = 80;
        
        NOISE_SEED = 135;
        with_noise = true;
        alphaAN = 0.1;
        noise;

        zc;
        zc_root = 25;
        zc_len  = 63;
        
        Nframe = 10;
        
        M = 4;
        phaseOffsetPSK = pi/4;

        rrc;
        data_symm;
        
    end
    
    methods
        function obj = QPSKTransmitter()
            obj.rrc = rcosdesign(obj.rolloff, obj.span, obj.os, 'sqrt');
            obj.zc = zadoffChuSeq(obj.zc_root, obj.zc_len);
        end
        

        function txSignal = transmit(obj)
            tx_sym = obj.buildFrameSymbols();

            tx_up = upsample(tx_sym, obj.os);
            tx_signal_frame = filter(obj.rrc, 1, tx_up);

            max_val = max(abs(tx_signal_frame));
            tx_signal_frame = 0.9 * (tx_signal_frame / max_val);
            % tx_signal_frame = (tx_signal_frame / max_val);
            txSignal = repmat(tx_signal_frame, obj.Nframe, 1);
        end
        

        function tx_sym = buildFrameSymbols(obj)
            % 1) Silence (4 symbols)
            silence1 = zeros(obj.ZEROS_LEN,1);
            
            % 2) Known constant symbols (+1)
            %known_sym = ones(obj.DUMMY_BIT_LEN,1);,
            known_bits = obj.generateSeededBits(obj.DUMMY_BIT_LEN, obj.M, obj.DUMMY_BIT_SEED);
            known_sym = obj.bits2psk(known_bits,obj.M, obj.phaseOffsetPSK);
            
            % 3) Zadoff - Chu Sequence
            zc_sym  = obj.zc;
            
            % 4) Plot symbols
            plot1_bits = obj.generateSeededBits(obj.PLOT_LEN, obj.M, obj.PLOT_SEED);
            plot1 = obj.bits2psk(plot1_bits,obj.M, obj.phaseOffsetPSK);
            obj.plot = plot1;

           switch lower(obj.data_mode)

                case 'visual'
                    % --------------------------------
                    % VISUAL MODE (TEXT BASED)
                    % --------------------------------
                    [data_bits, orgDataLen, ~] = obj.string2bits_linear('itü thal');
                    obj.ORG_DATA_LEN = orgDataLen;
                    obj.DATA_LEN = length(data_bits);
                    obj.DATA_SYMBOL_LEN = obj.DATA_LEN / log2(obj.M);
            
                case 'analysis'
                    % --------------------------------
                    % ANALYSIS MODE (SEEDED RANDOM)
                    % (COMMENT'Lİ KISIM AKTİF)
                    % --------------------------------
                    data_bits = obj.generateSeededBits(obj.DATA_LEN, obj.M, obj.DATA_SEED);
                    obj.ORG_DATA_LEN = length(data_bits);
                    obj.DATA_LEN = length(data_bits);
                    obj.DATA_SYMBOL_LEN = obj.DATA_LEN / log2(obj.M);
                otherwise
                    error('Unknown data_mode');
            end
            
            if obj.with_noise
                switch lower(obj.noise_mode)
                    case 'phase_shift_an'
                        data_sym = obj.bits2psk(data_bits, obj.M, obj.phaseOffsetPSK);
                        [noised, ~] = obj.addArtificialNoise(data_sym, obj.alphaAN, obj.NOISE_SEED);
                        data_sym = noised;
                    case 'gaussian_an'
                        data_sym = obj.bits2psk(data_bits, obj.M, obj.phaseOffsetPSK);
                        noised = obj.addGaussianArtificialNoise(data_sym, obj.alphaAN, obj.NOISE_SEED);
                        data_sym = noised;
                    case 'xor'
                        an = logical(obj.generateSeededBits(obj.DATA_LEN/2, obj.M, obj.NOISE_SEED));
                        data_bits_l = logical(data_bits);
                        noised = xor(data_bits_l, an);
                        noised = double(noised);   % sonra modülasyon için geri çevir
                        data_sym = obj.bits2psk(noised, obj.M, obj.phaseOffsetPSK);
                    otherwise
                        error('Unknown data_mode');
                end
            else
                data_sym = obj.bits2psk(data_bits, obj.M, obj.phaseOffsetPSK);
            end
            
            obj.data_symm = data_sym;
            obj.tx_data_bits = data_bits;

            % 6) Plot symbols
            %plot2 = sign(randn(obj.PLOT_LEN,1));
            plot2 = plot1;

            % 7) Silence
            silence2 = zeros(obj.ZEROS_LEN,1);
            
            % ---- FULL FRAME (SYMBOL LEVEL) ----
            tx_sym = [
                silence1;
                known_sym;
                zc_sym;
                plot1;
                data_sym;
                plot2;
                silence2
            ];
        end
    end

    methods (Static)
        
        function [y, phi] = addArtificialNoise(x, alphaAN, seed)
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
            y = x .* exp(1j * phi);
        end

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

        
        function [encodedBits, orgDataLen, G] = string2bits_linear(str)
            % STRING2BITS_LINEAR Convert string to bit sequence and apply Hamming (7,4) linear coding
            %
            % INPUT:
            %   str - input string
            % OUTPUT:
            %   encodedBits - linear coded bit sequence (column vector)
            %   G - generator matrix
        
            %% 1️⃣ String → bits
            ascii_vals = double(str);              % ASCII değerleri
            bits_mat = de2bi(ascii_vals, 8, 'left-msb'); % Her karakter 8-bit
            bits = bits_mat.';                     % sütun bazlı
            bits = bits(:);                        % column vector
            orgDataLen = length(bits);
            %% 2️⃣ Hamming (7,4) generator matrix
            G = [1 0 0 0 1 1 0;
                 0 1 0 0 1 0 1;
                 0 0 1 0 0 1 1;
                 0 0 0 1 1 1 1];
        
            %% 3️⃣ Bitleri 4'erli gruplara ayır ve padding
            n = 4;  % her grup 4 bit
            pad_len = mod(-length(bits), n);  % eksik bit sayısı
            if pad_len > 0
                bits = [bits; zeros(pad_len, 1)];
            end
            bits = reshape(bits, n, []);  % 4xN matris
        
            %% 4️⃣ Linear kodlama (4 bit → 7 bit)
            encodedBits = mod(G' * bits, 2);   % 7xN
            encodedBits = encodedBits(:);      % column vector
        end
        
        function y = addGaussianArtificialNoise(x, var, seed)
            
            rng(seed)
            len = length(x);
            an = sqrt(var/2) * (randn(len,1) + 1j*randn(len,1));
            y = x + an;

        end
        
    end
end