function j = getNollIndex(n,m)
    j1 = n*(n+1)/2;
    j2 = abs(m);    
    if mod(n,4) == 2 || mod(n,4) == 3    
        if m < 0
            j3 = 0;
        elseif m >= 0 
            j3 = 1;
        end
    elseif mod(n,4) == 0 || mod(n,4) == 1
        if m > 0
            j3 = 0;
        elseif m <= 0 
            j3 = 1;
        end
    end    
    j = j1 + j2 + j3 ;
end