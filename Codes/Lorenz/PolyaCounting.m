% Plot the number of UPOs for each sequence length 
% Written by: Prerna Patil 
% Date: 11th October 2023 

clc
clear 
close all 

LW = 'linewidth';
addpath('../utils')
% Hand Calculations  
Seqlength = 2:20; 
NumUPOs =   [1, 2, 3, 6, 9, 18, 30, 56, 99];

%% 
PolyaComputations = zeros(size(NumUPOs));
AllPossibleValues = 2.^Seqlength;
NumUPOsv2 = zeros(size(NumUPOs));
for i=1:length(Seqlength)
    N = Seqlength(i);
    PolyaComputations(i) = PolyaEnumerationTheorem(N,2);
    % Recounted 
    Recount = 2; % Subtract the monochromatic values 
    if(isprime(N))
        NumUPOsv2(i) = PolyaComputations(i) - Recount;
    else
        Divisors = 2:N-1;
        Divisors=Divisors(~(rem(N, Divisors)));
        for counter = 1:length(Divisors)
            Index = find(Seqlength==Divisors(counter));
            Recount = Recount + NumUPOsv2(Index);
        end
        NumUPOsv2(i) = PolyaComputations(i) - Recount;
    end
end

% Plot the values 
semilogy(Seqlength,NumUPOsv2,'*-k',LW,1.5)
hold on 
semilogy(Seqlength,PolyaComputations,'o-b',LW,1.5)
semilogy(Seqlength,AllPossibleValues,'d--r',LW,1.5)
xlabel('Symbolic Sequence length');
ylabel('Number of unique UPOs');
set(gca,'FontName','Times New Roman','FontSize',16)
legend('Unique UPOs','Polya Enumeration theorem','All possible combinations');

% Filename = sprintf('ViswanathFigures/NumberOfUniqueUPOs');
% saveas(gca,Filename,'epsc');
% saveas(gca,Filename,'fig');