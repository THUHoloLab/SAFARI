function [n,m] = getZernikeIndices(j)
    j_max = j;    
    n_max = floor(-0.5 + sqrt(0.25 + 2*j_max));    
    Z_ind = zeros(j_max,3);
    
    Z_ind(1,1) = 1;

    counter = 2;
    for n_step = 1:n_max
        for m_step = -n_step:2:n_step
            Z_ind(counter,1) = getNollIndex(n_step,m_step);
            Z_ind(counter,2) = n_step ;
            Z_ind(counter,3) = m_step ;        
            counter = counter + 1;
        end
    end    
    Z_ind = sortrows(Z_ind,1);    
    n = Z_ind(j,2);
    m = Z_ind(j,3);
end

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

        