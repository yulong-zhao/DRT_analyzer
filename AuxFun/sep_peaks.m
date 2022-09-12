function [x_sep] = sep_peaks(y, mode)
% This function seperates peaks for local attention.

[pks_i,locs_i,w_i,p_i] = findpeaks(y);

x_sep = [];


if (mode == 1)
    % seperate two peaks.
    
    x_begin = 1;
    
    
    while true
        locs_now = locs_i(locs_i>x_begin);
        
        if (length(locs_now)>1)
            x_end = locs_now(2);
            x_sep = [x_sep; x_begin, x_end];
            x_begin = locs_now(1);
        else
            x_sep = [x_sep; x_begin, length(y)];
            break;
        end
    end
    
elseif (mode == 2)
    % Try to isolate each peak.
    factor = 1.2;
    prev_peak = 1;
    next_peak = locs_i(2);
    for i=1:length(locs_i)
        
        if(i==1)
            x_begin = prev_peak;
            x_end   = min(next_peak, locs_i(i) + round(factor*w_i(i))  );
            
        elseif(i==length(locs_i))
            x_begin = max(prev_peak, locs_i(i) - round(factor*w_i(i))  );
            x_end   = next_peak;
            
        else
            x_begin = max(prev_peak, locs_i(i) - round(factor*w_i(i))  );
            x_end   = min(next_peak, locs_i(i) + round(factor*w_i(i))  );
        end
        
        
        

        
        x_sep = [x_sep; x_begin, x_end];
        
        prev_peak = locs_i(i);
        
        if( (i+1)>=length(locs_i) )
            next_peak = length(y);
        else
            next_peak = locs_i(i+2);
        end

    end
end




end