% function to generate an m-sequence via a LFSR
function seq = lfsr(seq_len, poly, init)

    % init lfsr w/ seed
    lfsr = init;

    for i = 1:seq_len
        seq(i) = lfsr(end);
        tmp = mod(sum(lfsr(poly)),2);
        lfsr(2:end) = lfsr(1:end-1);
        lfsr(1) = tmp;
    end
    
    seq = seq';
    
end
