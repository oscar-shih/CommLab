%% Start
clc;
clear all; %#ok<CLALL>
close all;
rng(0);
%% 2a encode data
bin = [1, 0, 1, 1, 0];
impulse_response = [1, 0, 0; 1, 0, 1; 1, 1, 1];
enc_data = conv_enc(bin, impulse_response);
fprintf('2a result: %s\n', strjoin(string(enc_data)));

%% 2b decode data
bin = [0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0];
impulse_response = [1, 0, 0; 1, 0, 1; 1, 1, 1];
dec_data = conv_dec(bin, impulse_response);
fprintf('2b result: %s\n', strjoin(string(dec_data)));

bin_1e = [0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0];
impulse_response_1e = [1, 0, 1; 1, 1, 1];
dec_data_1e = conv_dec(bin_1e, impulse_response_1e);
fprintf('1e result: %s\n', strjoin(string(dec_data_1e)));
%% 2c inference and plot
len = 100000;
bin = randi([0, 1], 1, len);
impulse_response = [1, 0, 0; 1, 0, 1; 1, 1, 1]; 
BER = zeros(1, 10);
ps = 0:0.1:1;
fprintf("2c result: \n");
for i = 1:length(ps)
    p = ps(i);
    enc_data = conv_enc(bin, impulse_response);
    Channel = channel(enc_data, p);
    dec_data = conv_dec(Channel, impulse_response);
    BER(i) = sum(dec_data ~= bin) / length(bin);
    fprintf('%f ', BER(i));
end
fprintf("\n");
fig_2c = plot(ps, BER);
title('BER - P');
xlabel("P");
ylabel("BER");
fig_2c(1).LineWidth = 2;
fig_2c(1).Color = 'green';
%% Problem3
len = 100000;
bin = randi([0, 1], 1, len);
impulse_response = [1, 1, 0; 1, 0, 1];
BER = zeros(1, 10);
ps = 0:0.1:1;
fprintf("Problem3 result: \n");
for i = 1:length(ps)
    p = ps(i);
    enc_data = conv_enc(bin, impulse_response);
    Channel = channel(enc_data, p);
    dec_data = conv_dec(Channel, impulse_response);
    BER(i) = sum(dec_data ~= bin) / length(bin);
    fprintf("%f ", BER(i));
end
fprintf("\n");
fig_3 = plot(ps, BER);
title('BER - P');
xlabel("P");
ylabel("BER");
fig_3(1).LineWidth = 2;
fig_3(1).Color = 'red'; 
%% Viterbi algorithm
function viterbi = create_viterbi(bin, impulse, table)
    len = length(bin) / size(impulse, 1);
    b_res = reshape(bin, size(impulse, 1), len).';
    viterbi = inf(len+1, 3, 4);
    viterbi(1, 2, 1) = 0;
    for j = 1:len
        for i = 1:4
            if viterbi(j, 2, i) ~= inf
                for k = 0:1
                    S = sum(xor(table(i+4*k, :), b_res(j, :)));
                    i_ = 2*k + floor((i-1)/2) + 1;
                    if viterbi(j, 2, i) + S < viterbi(j+1, 2, i_)
                        viterbi(j+1, :, i_) = [i,  viterbi(j, 2, i)+S,  k];
                    end
                end
            end
        end
    end
end
%% Construct table
function table = create_table(impulse)
    table = zeros(8, size(impulse, 1));
    for i = 0:7
        for j = 1:size(impulse, 1)
            S = sum([floor(i/4), floor(mod(i, 4)/2), mod(i, 2)].*impulse(j, :));
            table(i+1, j) = mod(S, 2);
        end
    end   
end
%% Encoding
function enc_data = conv_enc(bin, impulse)
    pre_state = 0;
    table = create_table(impulse);
    for i = 1:length(bin)
        cur_state = bin(i)*4 + pre_state;
        pre_state = floor(cur_state/2);
        S = size(impulse, 1) * (i - 1);
        enc_data(S + 1: S + size(impulse, 1)) = table(cur_state+1, 1:size(impulse, 1));
    end
end
%% Channel
function out = channel(enc, prob)
    mask = rand(1, length(enc)) < prob;
    out = xor(mask, enc);
end
%% Decoding
function dec_data = conv_dec(bin, impulse)
    cur_state = 1;
    table = create_table(impulse);
    viterbi = create_viterbi(bin, impulse, table);
    len = length(bin) / size(impulse, 1);
    for i = len: -1: 1
        dec_data(i) = viterbi(i+1, 3, cur_state);
        cur_state = viterbi(i+1, 1, cur_state);
    end
end