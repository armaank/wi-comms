% special bpsk decoder to map constellation to walsh-encoded values
function [demod_out] = bpskcdma(mod_in)

    demod_out = zeros(size(mod_in));
    
    for ii = 1:length(mod_in)
        
        if mod_in(ii) > 0.5
            demod_out(ii) = 1;
            
        elseif abs(mod_in(ii))< 0.5
            demod_out(ii) = 0;
            
        elseif mod_in(ii) < -0.5
            demod_out(ii) = -1;
            
        end
        
    end
end


