function [mL, g, H] = nestedfunU(v, lnp, Aeq, beq)
        K = size(v, 1);
        Aeq_ = Aeq';
        beq_= beq';
       
        lnx = lnp - 1 - Aeq_ * v;
        
        % robustificaton
        lnx = max(lnx, -150.0); 
        x = exp(lnx);   
        disp(x);
        % Lagrange dual function
        L = x' * (lnx - lnp + Aeq_ * v) - beq_ * v; 
        % take neg values since we want to maximize
        mL = -L;
        
        % gradient and Hessian
        g = beq - Aeq * x;    
        H = Aeq * ((x * ones(1, K)) .* Aeq_); 
end


