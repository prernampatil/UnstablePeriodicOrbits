function H = CreateHankelMatrix(xdat,height,width,skip)

H = zeros(height,width);
for k=1:height
    H(k,:) = xdat(k:skip:width*skip+k-1,1);
end
end