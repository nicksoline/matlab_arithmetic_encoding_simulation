%{
// ENCODING MARKOV BINARY SEQUENCE
//
// Nicholson Eugene
// Università degli studi di Siena
%}
fileID1 = fopen('Markov_sequence_Nicholson.bin');
eseq = fread(fileID1, 'ubit1');
fclose(fileID1);
eseq' % Print the sequence from the input file
Prob_calc = eseq';
eseq = eseq + 1; % adding 1 to all elements of the sequence to remove 0's
seq = eseq;
eseq = eseq';   % convert to row vector
A = unique(eseq); % sorting out repeated bits and storing only unique bits
h_enc = hist(eseq,A); % histogram of the bits
count = h_enc;

[row_seq, col_seq] = size(seq); 
if (row_seq > 1)
    seq = seq';
end

[row_count, col_count] = size(count); 
if (row_count > 1)
    count = count';
end

cum_count = [0, cumsum(count)];  % cumulative sum of histogram values
total_count = cum_count(end); 
N = ceil(log2(total_count)) + 2;  % round of the total bits to nearest integer
X = ceil(log2(length(count))) + 2;
LOW = 0; %initialize lower bound
HIGH = 2^N-1; %initialize upper bound
E3_count = 0; 
code_len = length(seq) * X + N;  % estimating length of code
code = zeros(1, code_len);   % initilize a row vector with zeros
code_index = 1;

%Calculating the probability of occurance of 1's and 0's
pc11 = 0;
pc10 = 0;
pc01 = 0;
pc00 = 0;
for pc = 1:length(seq)
    if pc > 999999 && Prob_calc(pc) == 0 && Prob_calc(pc-1) == 0
        pc00 = pc00 + 1;
    elseif pc > 999999 && Prob_calc(pc) == 0 && Prob_calc(pc-1) == 1
        pc01 = pc01 + 1;
    elseif pc > 999999 && Prob_calc(pc) == 1 && Prob_calc(pc-1) == 0
        pc10 = pc10 + 1;
    elseif pc > 999999 && Prob_calc(pc) == 1 && Prob_calc(pc-1) == 1
        pc11 = pc11 + 1;
    elseif Prob_calc(pc+1) == 1 && Prob_calc(pc) == 1
        pc11 = pc11 + 1;
    elseif Prob_calc(pc+1) == 1 && Prob_calc(pc) == 0
        pc10 = pc10 + 1;
    elseif Prob_calc(pc+1) == 0 && Prob_calc(pc) == 1
        pc01 = pc01 + 1;
    else
        pc00 = pc00 + 1;
    end
end
Prob00 = pc00/1000000;
Prob01 = pc01/1000000;
Prob10 = pc10/1000000;
Prob11 = pc11/1000000;
fprintf('The probability of a 0 from 0: %f\n',Prob00);
fprintf('The probability of a 0 from 1: %f\n',Prob01);
fprintf('The probability of a 1 from 0: %f\n',Prob10);
fprintf('The probability of a 1 from 1: %f\n\n',Prob11);

% encoding the sequence of bits
for k = 1:length(seq) 
    symbol = seq(k);   
    LOW_new = LOW + floor( (HIGH-LOW+1)*cum_count(symbol+1-1)/total_count ); % floor -> round towards negative infinity, cum_count -> cumulative sum
    HIGH = LOW + floor( (HIGH-LOW+1)*cum_count(symbol+1)/total_count )-1;
    LOW = LOW_new;
    
    while( isequal(bitget(LOW, N), bitget(HIGH, N)) || (isequal(bitget(LOW, N-1), 1) && isequal(bitget(HIGH, N-1), 0) ) )
        if isequal(bitget(LOW, N), bitget(HIGH, N))  
            bit = bitget(LOW, N);  %bitget(LOW, N) % returns bit value at position N to the integer array LOW
            code(code_index) = bit;
            code_index = code_index + 1;
            LOW = bitshift(LOW, 1) + 0;  % bits of LOW is shifted to the left by 1 bit & add 0 to LSB
            HIGH = bitshift(HIGH, 1) + 1; % bits of HIGH is shifted to the left by 1 bit & add 1 to LSB
%  %{           if (E3_count > 0)
%                 code(code_index:code_index + E3_count-1) = bitcmp(bit,1) * ones(1, E3_count);
%                 code_index = code_index + E3_count;
%                 E3_count = 0;
%             end
%            
            LOW = bitset(LOW, N+1, 0);  % N+1th bit of LOW is set to 0
            HIGH  = bitset(HIGH, N+1, 0);  % N+1th bit of HIGH is set to 0
            elseif ( (isequal(bitget(LOW, N-1), 1) && isequal(bitget(HIGH, N-1), 0) ) )
                LOW = bitshift(LOW, 1) + 0;   % 1st bit of LOW is set to 1
                HIGH  = bitshift(HIGH, 1) + 1;   % 1st bit of HIGH is set to 1
                LOW = bitset(LOW, N+1, 0);   % N+1thbit of LOW is set to 0
                HIGH = bitset(HIGH, N+1, 0);   % N+1th bit of HIGH is set to 0
                LOW = bitxor(LOW, 2^(N-1) );   % bit wise XOR of LOW with 2^(N-1)
                HIGH  = bitxor(HIGH, 2^(N-1) );   % bit wise XOR of HIGH with 2^(N-1)
%                E3_count = E3_count+1;
        end
    end
end

bin_LOW = de2bi(LOW, N, 'left-msb');
 if E3_count == 0
     code(code_index:code_index + N - 1) = bin_LOW;
     code_index = code_index + N;
 else
    bit = bin_LOW(1);
    code(code_index) = bit;
    code_index = code_index + 1;
    code(code_index:code_index+E3_count-1) = bitcmp(bit,1).*ones(1, E3_count);
    code_index = code_index + E3_count;
    code(code_index:code_index+N-2) = bin_LOW(2:N);
    code_index = code_index + N - 1;
end

code = code(1:code_index - 1);
code = code';

fileID6 = fopen('nicholson_encoder_code.bin','w');
fwrite(fileID6,code,'ubit1');
fclose(fileID6);

W = length(eseq);
X = length(count);
Z = length(code);

encoder_out = [W,X,Z];
fileID2 = fopen('nicholson_encodedseq.bin','w');
fwrite(fileID2,encoder_out(1),'ubit22');  % length of sequence
fwrite(fileID2,encoder_out(2),'ubit5');   % length of count
fwrite(fileID2,encoder_out(3),'ubit22');  % length of code
fwrite(fileID2,code,'ubit1');
fwrite(fileID2,count,'ubit22');
fclose(fileID2);
fprintf('\nENCODING COMPLETE\n');