function [Ns,r0_step]=  NumberOfScreens(rytov,k,z,criterion)
    % First find the C2n of the system
    c2n = rytov./( 1.23 * (k^(7/6)) * (z^(11/6)) );
    
    % Find the Fried parameter of the system
    r0 = (0.423 * c2n * k^2 * z)^(-3/5);
    r0_step = r0;
    
    disp(['The C2n number of the system is ',num2str(c2n)])
    disp(['Fried parameter is ',num2str(r0)])
    
    step_rytov = rytov;
    equiv_rytov = rytov;
    
    Ns = 1;
    
    while (step_rytov > criterion || equiv_rytov < rytov)
        Ns = Ns + 1;
        step_rytov = 5.32 * (k^(-5/6)) * (z^(5/6)) * (r0^(-5/3)) * ((1-1/(2*Ns))^(5/6)) / Ns;
        equiv_rytov = 0;
        for n = 1:Ns
            equiv_rytov = equiv_rytov + 5.32 * (k^(-5/6)) * (z^(5/6)) * (r0^(-5/3)) * ((1-(n-0.5)/Ns)^(5/6)) / Ns;
        end
    disp(["Rytov variance of first slab is ",num2str(step_rytov)])
    disp(["Number of screens required: ",num2str(Ns)])
    
    r0_step = (Ns^(3/5))*r0;
    end
end