function Combinations = PolyaEnumerationTheorem(N,m)
% m is number of colors 
% Number of symmetry groups 
G = N; 
Orbits = zeros(N,1);
% Compute the number of orbits for each symmtery group 
for i = 1:N
    Orbits(i) = gcd(N,i);
end
Combinations = 1/N*(sum(m.^Orbits));
