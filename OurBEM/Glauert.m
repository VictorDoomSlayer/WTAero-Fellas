function [aOut,Ct] = Glauert(aIn)
    Ct1 = 1.816;
    Ct2 = 2*sqrt(Ct1) - Ct1;
        if aIn > 0.5 % Glauert Correction for induction factors above 0.5
            Ct = 1.816 - 4*(sqrt(1.816)-1)*(1-aIn);
        else
            Ct = 4*aIn*(1-aIn);
        end

        if Ct >= Ct2
           aOut = 1 + (Ct - Ct1)/(4*sqrt(Ct1)-4);           
        else
           aOut = 0.5 - (sqrt(1-Ct)/2);
        end
end