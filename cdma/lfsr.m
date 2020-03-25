% function to generate an m-sequence via a LFSR
function seq = lfsr(seq_len, poly, init)

    % change form of polynomial to companion set
    order = poly(1);
    poly = [poly(1), fliplr(order - poly(2:end))];
    % init lfsr
    lfsr = init;
    seq = zeros(seq_len,1);
  
    % generate sequence
    for ii = 1:seq_len
        seq(ii) = lfsr(end);
        tmp = mod(sum(lfsr(poly)),2);
        lfsr(2:end) = lfsr(1:end-1);
        lfsr(1) = tmp;
    end
        
end
