function screen = NollPhaseScreen(N,D,r0,num_modes,rs,Zernikes)
   
    % Set MATLAB random number generator %
    rng(rs)
            
    % Generate the Noll Covariance Matrix %            
    NollCov = zeros(num_modes-1,num_modes-1);
    counter1 = 1;
    for j1 = 2:num_modes
        counter2 = 1;
        for j2 = 2:num_modes
            NollCov(counter1,counter2) = GenNollCovMatrix(j1,j2);    
            counter2 = counter2 + 1;
        end
        counter1 = counter1 + 1;
    end
    NollCov = NollCov * ((D/r0)^(5/3));    
    
    % Perform SVD and generate weighting vectors %        
    % Perform SVD %
    [~,~,Udag] = svd(NollCov);
    
    % Calculate weighting matrix %
    B = zeros(num_modes-1,1);    
    for B_ind = 1:num_modes-1
        var = NollCov(B_ind,B_ind);    
        sigma = sqrt(var);    
        B(B_ind) = sigma * randn(1);    
    end    
    A = Udag*B;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Superimpose Zernike Modes %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    screen = zeros(N,N);
    for j = 2:num_modes
        screen = screen + A(j-1) * Zernikes(:,:,j-1);
    end

end

function NollCov = GenNollCovMatrix(j1,j2)

    [n1,m1] = getZernikeIndices(j1);    
    [n2,m2] = getZernikeIndices(j2);

    if m1 == m2
        term1 = 2.2698 * sqrt((n1 + 1) * (n2 + 1)) * (-1)^((n1 + n2 - 2*m1)/2);
        term2 = gamma((n1 + n2 - 5/3)/2)/gamma((-n1+n2)/2+  17/6)/gamma(  (n1 - n2)/2 + 17/6    )/gamma( (n1+n2)/2  +23/6);
        NollCov = term1*term2;
    else
        NollCov = 0;
    end
end

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

