%{
// DECODING MARKOV BINARY SEQUENCE
//
// Nicholson Eugene
// Università degli studi di Siena
%}
fileID3 = fopen('nicholson_encodedseq.bin');
encoder_out(1) = fread(fileID3,1,'ubit22');
W = encoder_out(1); % length of sequence
encoder_out(2) = fread(fileID3,1,'ubit5');
X = encoder_out(2); % length of count
encoder_out(3) = fread(fileID3,1,'ubit22');
Z = encoder_out(3); % length of code
S = fread(fileID3,Z,'ubit1');
code = S';
h_dec = fread(fileID3,X,'ubit22');
count = h_dec';
len = W;
fprintf('\nENCODED FILE LOADED SUCCESFULLY\nPlease wait....\n');
[row_cd, col_cd] = size(code);
if (row_cd > 1)
    code = code';
end

[row_c, col_c] = size(count);
if (row_c > 1)
    count = count';
end

cum_counts = [0, cumsum(count)];  % cumulative sum of histogram values
total_count = cum_counts(end);   
N = ceil(log2(total_count)) + 2; % round of the total bits to nearest integer

low = 0;   
high = 2^N-1;
bin_tag = code(1:N);
dec_tag = bi2de(bin_tag, 'left-msb');
dec_seq = zeros(1,len);  % initilize a row vector of size 'len' with zeros
dec_seq_no = 1;
k = N; 
ptr = 0;

while (dec_seq_no <= len) 
    
    dec_tag_new = floor( ((dec_tag-low+1)*total_count-1)/(high-low+1) ); 
    % round to nearest integer less than or equal to the element, negative infinity (opp of ceil)
    
    if dec_tag_new == cum_counts(end)
    ptr = length(cum_counts)-1;
    return
    end
    
    c = find(cum_counts <= dec_tag_new);
    ptr = c(end);
    
    dec_seq(dec_seq_no) = ptr;  
    dec_seq_no = dec_seq_no + 1;
    low_new = low + floor( (high-low+1)*cum_counts(ptr-1+1)/total_count );  % floor -> round towards negative infinity, cum_count -> cumulative sum
    high = low + floor( (high-low+1)*cum_counts(ptr+1)/total_count ) - 1;
    low = low_new; 
    
    while ( isequal(bitget(low, N), bitget(high, N)) ||( isequal(bitget(low, N-1), 1) && isequal(bitget(high, N-1), 0) ) ),  
        
       if ( k == length(code) )
           break;
       end
       k = k + 1;
       if isequal(bitget(low, N), bitget(high, N))
            
            low = bitshift(low, 1) + 0; % bits of LOW is shifted to the left by 1 bit & add 0 to LSB
            high  = bitshift(high,  1) + 1;  % bits of HIGH is shifted to the left by 1 bit & add 1 to LSB
            dec_tag = bitshift(dec_tag, 1) + code(k);

            low = bitset(low, N+1, 0);  % N+1th bit of low is set to 0
            high  = bitset(high,  N+1, 0); % N+1th bit of HIGH is set to 0
            dec_tag = bitset(dec_tag, N+1, 0);
                  
        elseif ( isequal(bitget(low, N-1), 1) && isequal(bitget(high, N-1), 0) )
            low = bitshift(low, 1) + 0;   % 1st bit of LOW is set to 1
            high  = bitshift(high,  1) + 1; % 1st bit of HIGH is set to 1
            dec_tag = bitshift(dec_tag, 1) + code(k);
            
            
            low = bitset(low, N+1, 0);   % N+1th bit of low is set to 0
            high  = bitset(high,  N+1, 0); % N+1th bit of HIGH is set to 0
            dec_tag = bitset(dec_tag, N+1, 0);
            
            low = bitxor(low, 2^(N-1));  % bit wise XOR of LOW with 2^(N-1)
            high = bitxor(high,  2^(N-1)); % bit wise XOR of HIGH with 2^(N-1)
            dec_tag = bitxor(dec_tag, 2^(N-1));

        end
    end 
end

dec_seq = dec_seq-1 % Print the decoded Sequence
dec_seq = dec_seq';

fileID4 = fopen('nicholson_decoderseq.bin','w');
fwrite(fileID4,dec_seq,'ubit1');
fclose(fileID4);
fprintf('\nDECODING COMPLETE\n');