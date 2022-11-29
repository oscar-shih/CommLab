%% START b08502141
clc;
clear all; %#ok<CLALL>
close all;
rng(0);

%% 1a
symbol = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
prob = [0.2, 0.05, 0.005, 0.2, 0.3, 0.05, 0.045, 0.15];
H = 0;
for i=1:length(prob)
    H = H - prob(i)*log2(prob(i));
end
fprintf('H = %f\n', H);

result = huffman_dict(symbol, prob);
result  %#ok<NOPTS>
%% 1c
bits = {'10', '00010', '00000', '11', '01', '00001', '00001', '001'};
lengths = zeros(1, length(bits));
for i=1:length(bits)
    lengths(i) = length(bits{i});
end

ans_1c = 0;
for i=1:length(bits)
    ans_1c = ans_1c + 2^(-lengths(i));
end
fprintf('1c answer = %f\n', ans_1c);

%% 1d
L_bar = sum(prob.*lengths);
fprintf('L_bar = %f\n', L_bar);

%% 1e
target = 'gacab';
bitstream = '';
for i=1:length(target)
   for j=1:length(symbol)
       if target(i)==symbol{j}
           break
       end
   end
   bitstream = strcat( bitstream, bits{j});
end
fprintf('bitstream encoded: %s\n', bitstream);

%% 1f
decoded = '';
partial = '';
while ~isempty(bitstream)
    partial = strcat( partial, bitstream(1) );
    bitstream = bitstream(2:end);
    for i=1:length(bits)
        if strcmp(partial,bits{i})
            decoded = strcat( decoded, symbol{i} );
            partial = '';
            break
        end
    end
end
fprintf('test decoded: %s\n', decoded);

%% 1g
subset = cell(1, 15);
location = 1;
counter1 = 0;
counter2 = 0;
while counter1+counter2 < 1000 %% Random Sample for 1000 times
    candicate = rand(1,10);
    element = '';
    p = 1;
    for i=1:10
        symbol_index = 1;
        while sum(prob(1:symbol_index)) < candicate(i)
            symbol_index = symbol_index + 1;
        end
        element = strcat( element, symbol(symbol_index) );
        p = p * prob(symbol_index);
    end
    if abs(-log2(p)/10 - H) < 0.1
        counter1 = counter1 + 1;
        if counter1<=15
            subset{counter1} = element{1};
        end
    else
        counter2 = counter2 + 1;
    end
end
fprintf('probability: %f¢H\n', counter1/(counter1+counter2)*100);
for i=1:15
    fprintf('%s ', subset{i});
    if mod(i, 5)==0
        fprintf('\n');
    end
end

%% 2a -- as end of the file
symbols = { 's0', 's1', 's2', 's3', 's4' };
probs = [ 0.26, 0.25, 0.20, 0.15, 0.14 ];
example = huffman_dict(symbols, probs);
symbol = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
prob = [0.2, 0.05, 0.005, 0.2, 0.3, 0.05, 0.045, 0.15];
dict = huffman_dict(symbol, prob);
dict %#ok<NOPTS>

%% 2b
mydict = cell(15,5);
for i=1:8
    mydict{i,1} = symbol{i};
    mydict{i,2} = prob(i);
    mydict{i,5} = bits{i};
end
mydict{9,1} = 'cg'; mydict{9,2} = 0.05; mydict{9,3} = 3; mydict{9,4} = 7; mydict{9,5} = '0000';
mydict{10,1} = 'bf'; mydict{10,2} = 0.1; mydict{10,3} = 2; mydict{10,4} = 6; mydict{10,5} = '0001';
mydict{11,1} = 'cgbf'; mydict{11,2} = 0.15; mydict{11,3} = 9; mydict{11,4} = 10; mydict{11,5} = '000';
mydict{12,1} = 'cgbfh'; mydict{12,2} = 0.3; mydict{12,3} = 11; mydict{12,4} = 8; mydict{12,5} = '00';
mydict{13,1} = 'ad'; mydict{13,2} = 0.4; mydict{13,3} = 1; mydict{13,4} = 4; mydict{13,5} = '1';
mydict{14,1} = 'cgbfhe'; mydict{14,2} = 0.6; mydict{14,3} = 12; mydict{14,4} = 5; mydict{14,5} = '0';
mydict{15,1} = 'cgbfhead'; mydict{15,2} = 1.0; mydict{15,3} = 14; mydict{15,4} = 13;

encoded = huffman_enc({'g','a','c','a','b'}, mydict);
fprintf('encoded bits: %s\n', encoded);

%% 2c

decoded = huffman_dec(encoded, mydict);
fprintf('decoded bits: ');
disp(decoded);

%% 3a
candicate = rand(1,10);
symbols_selected = cell(1,10);
for i=1:10
    index = 1;
    while sum(prob(1:index)) < candicate(i)
        index = index + 1;
    end
    symbols_selected{i} = symbol{index};
end
binary = huffman_enc(symbols_selected, mydict);

fprintf('symbol sequence: %s\n', [symbols_selected{:}]);
fprintf('binary data: %s\n', binary);
fprintf('binary data length: %d\n', length(binary));

%% 3b
n = 10;
R = 200;

candicate = rand(R, n);
bin_length = zeros(1, R);
for k=1:R
    symbols_selected = cell(1, n);
    for i=1:n
        index = 1;
        while sum(prob(1:index)) < candicate(k, i)
            index = index + 1;
        end
        symbols_selected{i} = symbol{index};
    end
    binary = huffman_enc(symbols_selected, mydict);
    bin_length(k) = length(binary);
end
length_mean = mean(bin_length);

figure;
histogram(bin_length);
title(sprintf('$L_{%d}$(%d) = %f', n, R, length_mean), 'interpreter', 'latex');
xlabel('binary length (bits)');
ylabel('occurence (times)');

%% 3c
Ns = [10, 50, 100];
curves = cell(1, 3);

for N_index=1:length(Ns)
    n = Ns(N_index);
    Rs = [10, 20, 50, 100, 200, 500, 1000];
    L_bars = zeros(1, length(Rs));
    for R_index=1:length(Rs)
        R = Rs(R_index);
        candicate = rand(R, n);
        bin_length = zeros(1, R);
        for k=1:R
            symbols_selected = cell(1, n);
            for i=1:n
                index = 1;
                while sum(prob(1:index)) < candicate(k, i)
                    index = index + 1;
                end
                symbols_selected{i} = symbol{index};
            end
            binary = huffman_enc(symbols_selected, mydict);
            bin_length(k) = length(binary);
        end
        L_bars(R_index) = mean(bin_length)/n;
    end
    curves{N_index} = L_bars;
end

figure;
semilogx(Rs, curves{1}, 'red', Rs, curves{2}, 'green', Rs, curves{3}, 'blue');
text(Rs(end)+5, curves{1}(end), 'n = 10', 'color', 'red');
text(Rs(end)+5, curves{2}(end), 'n = 50', 'color', 'green');
text(Rs(end)+5, curves{3}(end)-0.008, 'n = 100', 'color', 'blue');


yL1 = yline(L_bar, '--', sprintf('$\\bar{L}$ = %f', L_bar), 'interpreter', 'latex');
yline(H, '--', sprintf('\\textsl{H}[\\textsl{X}] = %f', H), 'interpreter', 'latex');
xlabel('R (times)');
ylabel('$\bar{L}$', 'interpreter', 'latex');
title('$\bar{L} - R$', 'interpreter', 'latex');

%% BONUS
fprintf('\nBOUNS: ---------------\n');
symbol_old = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
prob_old = [0.2, 0.05, 0.005, 0.2, 0.3, 0.05, 0.045, 0.15];
symbol = cell(1, 8^3);
prob = zeros(1, 8^3);
index = 1;
for i=1:8
    for j=1:8
        for k=1:8
            symbol{index} = strcat( symbol_old{i}, symbol_old{j}, symbol_old{k} );
            prob(index) = prob_old(i) * prob_old(j) * prob_old(k);
            index = index + 1;
        end
    end
end

H = 0;
for i=1:length(prob)
    H = H - prob(i)*log2(prob(i));
end
H = H/3;
fprintf('H = %f\n', H);

%% bonus - 2a
dict = huffman_dict(symbol, prob);
% dict(1:512, :) %#ok<NOPTS>

L_bar = 0;
index = 1;
while isempty(dict{index, 3})
    L_bar = L_bar + dict{index, 2} * length(dict{index, 5}) / 3;
    index = index + 1;
end
fprintf('L_bar = %f\n', L_bar);

%% bonus - 2b

encoded = huffman_enc({'gac','aba'}, dict);
fprintf('encoded bits: %s\n', encoded);

%% bonus - 2c

decoded = huffman_dec(encoded, dict);
fprintf('decoded bits: ');
disp(decoded);

%% bonus - 3a
candicate = rand(1,10);
symbols_selected = cell(1,10);
for i=1:10
    index = 1;
    while sum(prob(1:index)) < candicate(i)
        index = index + 1;
    end
    symbols_selected{i} = symbol{index};
end
binary = huffman_enc(symbols_selected, dict);

fprintf('symbol sequence: %s\n', [symbols_selected{:}]);
fprintf('binary data: %s\n', binary);
fprintf('binary data length: %d\n', length(binary));

%% bonus - 3b
n = 10;
R = 200;

candicate = rand(R, n);
bin_length = zeros(1, R);
for k=1:R
    symbols_selected = cell(1, n);
    for i=1:n
        index = 1;
        while sum(prob(1:index)) < candicate(k, i)
            index = index + 1;
        end
        symbols_selected{i} = symbol{index};
    end
    binary = huffman_enc(symbols_selected, dict);
    bin_length(k) = length(binary);
end
length_mean = mean(bin_length)/3;

figure;
histogram(bin_length);
title(sprintf('$L_{%d}$(%d) = %f', n, R, length_mean), 'interpreter', 'latex');
xlabel('binary length (bits)');
ylabel('occurence (times)');

%% bonus - 3c
Ns = [10, 50, 100];
curves = cell(1, 3);

for N_index=1:length(Ns)
    n = Ns(N_index);
    Rs = [10, 20, 50, 100, 200, 500, 1000];
    L_bars = zeros(1, length(Rs));
    for R_index=1:length(Rs)
        R = Rs(R_index);
        candicate = rand(R, n);
        bin_length = zeros(1, R);
        for k=1:R
            symbols_selected = cell(1, n);
            for i=1:n
                index = 1;
                while sum(prob(1:index)) < candicate(k, i)
                    index = index + 1;
                end
                symbols_selected{i} = symbol{index};
            end
            binary = huffman_enc(symbols_selected, dict);
            bin_length(k) = length(binary);
        end
        L_bars(R_index) = mean(bin_length)/n/3;
    end
    curves{N_index} = L_bars;
end

figure;
semilogx(Rs, curves{1}, 'red', Rs, curves{2}, 'green', Rs, curves{3}, 'blue');
text(Rs(end)+5, curves{1}(end), 'n = 10', 'color', 'red');
text(Rs(end)+5, curves{2}(end), 'n = 50', 'color', 'green');
text(Rs(end)+5, curves{3}(end), 'n = 100', 'color', 'blue');


yL1 = yline(L_bar, '--', sprintf('$\\bar{L}$ = %f', L_bar), 'interpreter', 'latex');
yline(H, '--', sprintf('\\textsl{H}[\\textsl{X}] = %f', H), 'interpreter', 'latex');
xlabel('R (times)');
ylabel('$\bar{L}$', 'interpreter', 'latex');
title('$\bar{L} - R$', 'interpreter', 'latex');

%% FUNCTIONS ==========================================================================
%% huffman_dict function
function result = huffman_dict(symbols, prob)
    n = length(symbols);

    queue = cell(n, 3);
    for i=1:n
        queue{i,1} = prob(i);
        queue{i,2} = symbols{i};
        queue{i,3} = i;
    end
      
    indices = 1:n;
    result = cell(2*n-1, 5);    
    
    for i=1:n
       indices(i) = i;
       result{i, 1} = symbols{i};
       result{i, 2} = prob(i);
       result{i, 5} = '';
    end
    
    location = n+1;
    while size(queue,1)>1
        queue = sortrows(queue, [1,2,3]);
        temp = sortrows(queue([1,2], :), [2,1,3]);
        
        if strcmp(temp{1,2}, queue{2,2})
            queue([1,2], :) = queue([2,1], :);
        end
        result{location, 1} = strcat( result{queue{1,3}, 1}, result{queue{2,3}, 1} );
        result{location, 2} = result{queue{1,3}, 2} + result{queue{2,3}, 2};
        result{location, 3} = queue{1,3};
        result{location, 4} = queue{2,3};
        result{location, 5} = '';
      
        thesize = size(queue, 1);
        queue{thesize+1, 1} = result{location, 2}; 
        queue{thesize+1, 2} = result{location, 1}; 
        queue{thesize+1, 3} = location;            
        queue = queue(3:end, 1:end);
        location = location + 1;
    end
    
    for i=(2*n-1):-1:(n+1)
        L_child = result{i, 3};
        R_child = result{i, 4};
        result{L_child, 5} = strcat(result{i, 5}, '0');
        result{R_child, 5} = strcat(result{i, 5}, '1');
    end
end

%% huffman_enc function
function bin_seq = huffman_enc(sym_seq, dict)
    bin_seq = '';
    
    boundary = 1;
    while isempty(dict{boundary,3})
        boundary = boundary + 1;
    end
    boundary = boundary - 1;
    
    for i=1:length(sym_seq)
        for j=1:boundary
            if strcmp(dict{j,1}, sym_seq{i})
                bin_seq = strcat( bin_seq, dict{j,5} );
                break
            end
        end
    end
end

%% huffman_dec function
function sym_seq = huffman_dec(bin_seq, dict)
    sym_seq = {};
    partial = '';
    
    boundary = 1;
    while isempty(dict{boundary,3})
        boundary = boundary + 1;
    end
    boundary = boundary - 1;
    
    while ~isempty(bin_seq)
        partial = strcat( partial, bin_seq(1) );
        bin_seq = bin_seq(2:end);
        for i=1:boundary
            if strcmp(partial, dict{i, 5})
                sym_seq{end+1} = dict{i, 1}; %#ok<AGROW>
                partial = '';
                break
            end
        end
    end
end