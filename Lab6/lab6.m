%% START
clc;
clear all; %#ok<CLALL>
close all;
rng(0);

%% 1
test({round(rand(1,15)),  2, 2, 'PAM'});
test({round(rand(1,15)),  4, 2, 'PAM'});
test({round(rand(1,15)),  8, 2, 'PAM'});
test({round(rand(1,15)),  16, 2, 'PAM'});

test({round(rand(1,10)),  2, 3, 'PSK'});
test({round(rand(1,10)),  4, 3, 'PSK'});
test({round(rand(1,10)),  8, 3, 'PSK'});
test({round(rand(1,10)),  16, 3, 'PSK'});
% 
test({round(rand(1,20)), 4, 4, 'QAM'});
test({round(rand(1,20)), 16, 4, 'QAM'});
test({round(rand(1,20)), 64, 4, 'QAM'});

%% 2a
d_min = 1;
M = 4;
bit_length = 10^4;
b = round(rand(1, bit_length));
u = symbol_mapper(b, M, d_min, 'PSK');
figure;
histogram2(u(:,1), u(:,2), 100);
hold on;
line([0, 0], [-0.5,0.5], 'LineWidth', 2, 'Color', 'r');
line([-0.5, 0.5], [0,0], 'LineWidth', 2, 'Color', 'r');

Eb = (d_min/2)^2 / log2(M) / (sin(pi/M))^2;
for N0 = [Eb/1, Eb/10,  Eb/100]
    n = normrnd(0 , N0, [bit_length/log2(M), 2]);
    r = u + n;
    figure;
    h = histogram2(r(:,1), r(:,2), 100);
    hold on;
    x_range = max(abs(min(r(:, 1))), abs(max(r(:, 1))));
    y_range = max(abs(min(r(:, 2))), abs(max(r(:, 2))));
    line([0, 0], [-y_range, y_range], 'LineWidth', 2, 'Color', 'r');
    line([-x_range, x_range], [0, 0], 'LineWidth', 2, 'Color', 'r');
    
    if N0==Eb/1     % keep r
        r1 = r;
    elseif N0==Eb/10
        r2 = r;
    else
        r3 = r;
    end
        
end

%% 2b
u0 = symbol_decider(u, M, d_min, 'PSK');
u1 = symbol_decider(r1, M, d_min, 'PSK');
u2 = symbol_decider(r2, M, d_min, 'PSK');
u3 = symbol_decider(r3, M, d_min, 'PSK');
fprintf('2b result: symbol error rate = %d, %d, %d, %d (respectively)\n', mean(sum(u'==u0')~=size(u, 2)), mean(sum(u'==u1')~=size(u, 2)), mean(sum(u'==u2')~=size(u, 2)), mean(sum(u'==u3')~=size(u, 2)));

%% 3a
name = 'PAM';
test_x = 0:0.5:10;
real_y = zeros(1, length(test_x));
figure;
color = {'red', 'green', 'blue', 'black'};
for M = [2,4,8,16]
    SNRs = 0:10;    % dB
    SERs = zeros(1, 11);
    for i=1:11 
        SNR = SNRs(i);
        d_min = 1;
        Eb = d_min^2 / 12 / log2(M) * (M^2 - 1); 
        N0 = Eb / 10^(SNR/10);
        bit_length = log2(M) * 10^7;    
        
        b = round(rand(1, bit_length));
        u = symbol_mapper(b, M, d_min, name);
        n = normrnd(0 , sqrt(N0/2), [bit_length/log2(M), 1]); 
        r = u + n;
        u_decided = symbol_decider(r, M, d_min, name);
        SERs(i) = mean(u'~=u_decided');                   
    end
    
    for i=1:length(real_y)
        real_y(i) = 2*qfunc(sqrt(6*log2(M)/(M^2-1) * (10^(test_x(i)/10)))); 
    end
    semilogy(SNRs, SERs, 'o', 'color', color{log2(M)});
    hold on;
    semilogy(test_x, real_y, '-', 'color', color{log2(M)});
    hold on;
end
title('SER - SNR for PAM');
xlabel('SNR (dB)');
ylabel('SER');
legend('M=2 Simulation', 'M=2 Theoretical', 'M=4 Simulation', 'M=4 Theoretical', 'M=8 Simulation', 'M=8 Theoretical', 'M=16 Simulation', 'M=16 Theoretical');
hold off;
fprintf('3a finished\n');

%% 3b
name = 'PSK';
test_x = 0:0.5:10;
real_y = zeros(1, length(test_x));
figure;
color = {'red', 'green', 'blue', 'black'};
for M = [2,4,8,16]
    SNRs = 0:10; 
    SERs = zeros(1, 11);
    for i=1:11 
        SNR = SNRs(i);
        d_min = 1;
        Eb = (d_min/2)^2 / log2(M) / (sin(pi/M))^2;
        N0 = Eb / 10^(SNR/10);
        bit_length = log2(M) * 10^7; 
        
        b = round(rand(1, bit_length));
        u = symbol_mapper(b, M, d_min, name);
        n = normrnd(0, sqrt(N0/2), [bit_length/log2(M), 2]);       
        r = u + n;
        u_decided = symbol_decider(r, M, d_min, name);
        SERs(i) = mean(sum(u'==u_decided')~=size(u, 2));    
    end
    
    for i=1:length(real_y)
        real_y(i) = 2*qfunc(sqrt(2*log2(M)*((sin(pi/M))^2) * (10^(test_x(i)/10))));
    end
    semilogy(SNRs, SERs, 'o', 'color', color{log2(M)});
    hold on;
    semilogy(test_x, real_y, '-', 'color', color{log2(M)});
    hold on;
end
title('SER - SNR for PSK');
xlabel('SNR (dB)');
ylabel('SER');
legend('M=2 Simulation', 'M=2 Theoretical', 'M=4 Simulation', 'M=4 Theoretical', 'M=8 Simulation', 'M=8 Theoretical', 'M=16 Simulation', 'M=16 Theoretical');
hold off;
fprintf('3b finished\n');

%% 3c
name = 'QAM';
test_x = 0:0.5:10;
real_y = zeros(1, length(test_x));
figure;
color = {'red', 'green', 'blue'};
for M = [4,16,64]
    SNRs = 0:10; 
    SERs = zeros(1, 11);
    for i=1:11 
        SNR = SNRs(i);
        d_min = 1;
        Eb = d_min^2 / 6 / log2(M) * (M - 1);
        N0 = Eb / 10^(SNR/10);
        bit_length = log2(M) * 10^7; 
        
        b = round(rand(1, bit_length));
        u = symbol_mapper(b, M, d_min, name);
        n = normrnd(0, sqrt(N0/2), [bit_length/log2(M), 2]);       
        r = u + n;
        u_decided = symbol_decider(r, M, d_min, name);
        SERs(i) = mean(sum(u'==u_decided')~=size(u, 2));    
    end
    
    for i=1:length(real_y)
        real_y(i) = 4*qfunc(sqrt(3*log2(M)/(M - 1) * (10^(test_x(i)/10))));
    end
    semilogy(SNRs, SERs, 'o', 'color', color{log2(M)/2});
    hold on;
    semilogy(test_x, real_y, '-', 'color', color{log2(M)/2});
    hold on;
end
title('SER - SNR for QAM');
xlabel('SNR (dB)');
ylabel('SER');
legend('M=4 Simulation', 'M=4 Theoretical', 'M=16 Simulation', 'M=16 Theoretical', 'M=64 Simulation', 'M=64 Theoretical');
hold off;
fprintf('3c finished\n');

%% 4a
name = 'PSK';
M = 2;
impulse_response = [1, 0, 1; 1, 1, 1];
figure;

SNRs = 0:0.5:5;
BERs = zeros(1, 11);
for i=1:11 
    SNR = SNRs(i);
    d_min = 1;
    Eb = (d_min/2)^2;
    N0 = Eb / 10^(SNR/10);
    bit_length = log2(M) * 10^6; 
    
    b = round(rand(1, bit_length));
    encoded_data = conv_enc(b, impulse_response);
    u = symbol_mapper(encoded_data, M, d_min, name);
    
    n = normrnd(0, sqrt(N0/2), [size(u, 1), size(u, 2)]);       
    r = u + n;
    
    u_decoded = MD_symbol_demapper(r, M, d_min, name);
    decoded_data = conv_dec(u_decoded, impulse_response);
    
    BERs(i) = mean(b~=decoded_data);    
end
semilogy(SNRs, BERs, '-o', 'color', 'red');
hold on;
fprintf('4a finished\n');

%% 4b
name = 'PSK';
M = 2;
impulse_response = [1, 0, 1; 1, 1, 1];

SNRs = 0:0.5:5;    % dB
BERs = zeros(1, 11);
for i=1:11 
    SNR = SNRs(i);
    d_min = 1;
    Eb = (d_min/2)^2;
    N0 = Eb / 10^(SNR/10);
    bit_length = log2(M) * 10^6;
    
    b = round(rand(1, bit_length));
    encoded_data = conv_enc(b, impulse_response);
    u = symbol_mapper(encoded_data, M, d_min, name);
    
    n = normrnd(0, sqrt(N0/2), [size(u, 1), size(u, 2)]);       
    r = u + n;
    
    decoded_data = soft_dec(r, impulse_response, M, d_min, name);
    
    BERs(i) = mean(b~=decoded_data);    
end
semilogy(SNRs, BERs, '-x', 'color', 'green');

title('BER - SNR for BPSK');
xlabel('SNR (dB)');
ylabel('BER');
legend('hard decoding', 'soft decoding');
hold off;
fprintf('4b finished\n');

%% FUNCTIONS
function mapping = get_mapping(M, d, name)
    grayCodes = {{'0' '1'}, {'00' '01' '11' '10'}, {'000' '001' '011' '010' '110' '111' '101' '100'}, {'0000' '0001' '0011' '0010' '0110' '0111' '0101' '0100' '1100' '1101' '1111' '1110' '1010' '1011' '1001' '1000'}};
    if strcmp(name, 'PAM')
        Ep = (d/2)^2;
        graycode = bin2dec(grayCodes{log2(M)});
        mapping = zeros(M, 1);
        for m=1:M
            mapping( graycode(m) + 1, 1) = (2*m - 1 - M) * (Ep^0.5);
        end
    elseif strcmp(name, 'PSK')
        Ep = (d/2/sin(pi/M))^2;
        graycode = bin2dec(grayCodes{log2(M)});
        mapping = zeros(M, 2);
        if M == 4
            const = 3*pi/4;
        else
            const = 0;
        end
        
        for m=1:M
            mapping( graycode(m) + 1, 1) = cos(2*pi*(m-1)/M + const) * (Ep^0.5);
            mapping( graycode(m) + 1, 2) = sin(2*pi*(m-1)/M + const) * (Ep^0.5);
        end
    elseif strcmp(name, 'QAM')
        Ep = (d/2)^2;
        graycode = bin2dec(grayCodes{log2(M^0.5)});
        mapping = zeros(M, 2);
        
        for head = 1 : (M^0.5)
            for tail = 1 : (M^0.5)
                key = graycode(head) * (M^0.5) + graycode(tail) + 1;
                mapping(key, 1) = (2*head - 1 - M^0.5) * (Ep^0.5);
                mapping(key, 2) = (2*tail - 1 - M^0.5) * (Ep^0.5);
            end
        end        
    else
        fprintf("Error: the name of the modulation should be 'PAM', 'PSK' or 'QAM'\n");
    end
end

%% symbol mapper
function sym_seq = symbol_mapper(bin_seq, M, d, name)
    mapping = get_mapping(M, d, name);
    
    if mod(length(bin_seq), log2(M))~=0
       fprintf("Warning: The length of the binary sequence should be multiples of %d. 0's are padded at the end.\n", log2(M));
       bin_seq(end+1: ceil(length(bin_seq)/log2(M))*log2(M)) = 0;
    end
    
    if strcmp(name, 'PAM')
        sym_seq = zeros(length(bin_seq)/log2(M), 1);
    else
        sym_seq = zeros(length(bin_seq)/log2(M), 2);
    end
    
    for i=1:length(bin_seq)/log2(M)
        key = 0;
        for j=1:log2(M)
            key = key*2 + bin_seq( (i-1)*log2(M) + j );
        end
        sym_seq(i, :) = mapping(key+1, :);
    end
end

%% test symbol mapper
function result = test(testcase)
    fprintf('test data: %s, M = %d, d = %d, module = %s\nresult: ', strjoin(string(testcase{1}), ''), testcase{2}, testcase{3}, testcase{4});
    
    result = symbol_mapper(testcase{1}, testcase{2}, testcase{3}, testcase{4});
    for i=1:length(result)
        if size(result, 2)==1
            fprintf('%f ', result(i, 1));
        else
            fprintf('(%f, %f) ', result(i, 1), result(i, 2));
        end
    end
    fprintf('\n\n');
end

%% MD symbol demapper
function bin_seq = MD_symbol_demapper(sym_seq, M, d, name)
    mapping = get_mapping(M, d, name);
    
    dec_seq = dsearchn(mapping, sym_seq) - 1;
    bin_seq = dec2bin(dec_seq, log2(M));
    bin_seq = bin_seq=='1';
    bin_seq = reshape( bin_seq' , 1 , size(bin_seq, 1)*size(bin_seq, 2));
end

%% test SER
function symbols = symbol_decider(sym_seq, M, d, name)
    mapping = get_mapping(M, d, name);
    
    result = dsearchn(mapping, sym_seq);
    symbols = mapping(result, :);
end

%% build finite state machine
function [postStates,Outputs] = build_FSM(impulse_response)
    n = size(impulse_response, 1);
    K = size(impulse_response, 2);
    
    postStates = zeros(1, 2^K);
    Outputs = zeros(2^K, n);
    states = dec2bin(0: 2^(K-1)-1);
    for i=1:2^(K-1) 
        state = states(i, :);
        for input=0:1
            status = zeros(1, K);
            status(1) = input;
            for j=1:K-1
                if state(j)=='1'
                    status(j+1) = 1;
                end
            end
            output = zeros(1, n);
            for j=1:n
                output(j) = mod(sum(    status .* impulse_response(j, :)     ),2);
            end
            poststate = bin2dec(strcat( string(input), state(1:end-1) )) + 1;
            
            key = bin2dec(strcat(state, string(input))) + 1;
            postStates(key) = poststate;
            Outputs(key, :) = output;
        end
    end
end

%% convolutional encoding
function encoded_data = conv_enc(binary_data, impulse_response)
    n = size(impulse_response, 1);   
    [postStates,Outputs] = build_FSM(impulse_response);
    
    state = 1;
    encoded_data = zeros(1, n*length(binary_data));
    for i=1:length(binary_data)
        key = (state-1)*2 + binary_data(i) + 1;
        state = postStates(key);   
        output = Outputs(key, :);
        
        for j=1:n
            encoded_data( (i-1)*n + j ) = output(j);
        end
    end
end

%% convolutional decoding
function decoded_data = conv_dec(binary_data, impulse_response)
    n = size(impulse_response, 1);
    K = size(impulse_response, 2);    
    [postStates,Outputs] = build_FSM(impulse_response);
    
    if mod(length(binary_data), n)~=0
       fprintf("Warning: The length of the binary data should be multiples of %d. 0's are padded at the end.\n", n);
       binary_data(end+1: ceil(length(binary_data)/n)*n) = 0;
    end

    T = length(binary_data)/n;
    output_seq = reshape(binary_data, n, T)';
    
    INF = length(binary_data) * 10 + 100;
    
    Dtable = ones(T+1, 2^(K-1))*INF;           
    Ptable = zeros(T+1, 2^(K-1));               
    Itable = zeros(T+1, 2^(K-1)); 
    isPossible = zeros(1, 2^(K-1));
    
    Dtable(1, 1) = 0;
    isPossible(1) = 1;
    
    for i=1:T
        isPossible_new = zeros(1, 2^(K-1));
        for j=1:2^(K-1)
            if ~isPossible(j) 
                continue;
            end
            
            for input=0:1
                key = (j-1)*2 + input + 1;
                nextstate = postStates(key);
                output = Outputs(key, :);
                HammingDistance = sum(output ~= output_seq(i, :));
                
                if HammingDistance + Dtable(i, j) < Dtable(i+1, nextstate)
                    Dtable(i+1, nextstate) = HammingDistance + Dtable(i, j);
                    Ptable(i+1, nextstate) = j;
                    Itable(i+1, nextstate) = input;
                end
                isPossible_new(nextstate) = 1;
            end
        end
        isPossible = isPossible_new;
    end
    
    % trace back
    miniState = 0;
    miniValue = INF;
    for i=1:2^(K-1)
        if Dtable(end, i) < miniValue
            miniValue = Dtable(end, i);
            miniState = i;
        end
    end
    
    stateSelected = miniState;
    decoded_data = zeros(1, T);
    for i = T : -1 : 1
        decoded_data(i) = Itable(i+1, stateSelected);
        stateSelected = Ptable(i+1, stateSelected);
    end
end

%% soft decoding
function decoded_data = soft_dec(symbol_data, impulse_response, M, d, name)
    mapping = get_mapping(M, d, name);
    [postStates,Outputs] = build_FSM(impulse_response);
    
    n = size(impulse_response, 1);
    K = size(impulse_response, 2);
    T = size(symbol_data, 1) * log2(M) / n;    
    
    Dtable = ones(T+1, 2^(K-1))*inf;
    Ptable = zeros(T+1, 2^(K-1));  
    Itable = zeros(T+1, 2^(K-1));       
    isPossible = zeros(1, 2^(K-1));
    
    Dtable(1, 1) = 0;
    isPossible(1) = 1;
    
    for i=1:T
        isPossible_new = zeros(1, 2^(K-1));
        for j=1:2^(K-1)
            if ~isPossible(j) 
                continue;
            end
            
            for input=0:1
                key = (j-1)*2 + input + 1;
                nextstate = postStates(key);
                output_c = Outputs(key, :);
                output_c = reshape(output_c, log2(M), n/log2(M))';
                EuclideanDistance = 0;
                for k=1:length(output_c)/log2(M)
                    point1 = mapping(bin2dec(num2str(output_c(k, :))) + 1, :);
                    point2 = symbol_data((i-1)*n/log2(M)+k, :);
                    EuclideanDistance = EuclideanDistance + norm(point1 - point2);
                end
                
                if EuclideanDistance + Dtable(i, j) < Dtable(i+1, nextstate)
                    Dtable(i+1, nextstate) = EuclideanDistance + Dtable(i, j);
                    Ptable(i+1, nextstate) = j;
                    Itable(i+1, nextstate) = input;
                end
                isPossible_new(nextstate) = 1;
            end
        end
        isPossible = isPossible_new;
    end
    
    miniState = 0;
    miniValue = inf;
    for i=1:2^(K-1)
        if Dtable(end, i) < miniValue
            miniValue = Dtable(end, i);
            miniState = i;
        end
    end
    
    stateSelected = miniState;
    decoded_data = zeros(1, T);
    for i = T : -1 : 1
        decoded_data(i) = Itable(i+1, stateSelected);
        stateSelected = Ptable(i+1, stateSelected);
    end
end