function [cycle,pass] = compute_Cycle_Pass(abs_orbit,asc_flag,lat_start)

    switch asc_flag % A means ascending, D descending
        case 'A'
            if lat_start<0 % ascending south of equator
                P = +1;
            else % ascending north of equator
                P = -1;
            end
        case 'D'
            P = 0;
    end
    P = P + abs_orbit*2 - 19;
    Q = mod(P,10688);
    R = mod(Q,2462);
    
    cycle = floor(P/10688)*13 + floor(Q/2462)*3 + floor(R/840) + 1;  % resulting cycle number
    pass = mod(R,840) + 1;                   	% resulting pass number

    
end

