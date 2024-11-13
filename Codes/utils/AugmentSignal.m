function SignalAugmented = AugmentSignal(Signal,Period,TotalLength)
nVars = size(Signal,2);
SignalAugmented = zeros(TotalLength,nVars);

for counter= 1:TotalLength
    Index = rem(counter,Period);
    if(Index==0)
        Value = Signal(end,:);
    else
        Value = Signal(Index,:);
    end
    SignalAugmented(counter,:) = Value;
end
