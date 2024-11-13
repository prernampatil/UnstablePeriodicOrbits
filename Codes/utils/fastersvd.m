function [U0, S0, V0] = fastersvd(A)

R1 = size(A,1); 
R2 = size(A,2); 

if(R1<R2)
    A2 = A*A'; 
    [U0,S2,~] = svd(A2); 
    S0 = diag(sqrt(diag(S2)));
    V0 = diag(1./diag(S0))*U0'*A; 
    V0 = V0';
    
else
    A2 = A'*A; 
    [V0,S2,~] = svd(A2);
    S0 = diag(sqrt(diag(S2)));
    U0 = A*V0*diag(1./diag(S0));
end