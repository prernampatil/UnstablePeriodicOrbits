 function s = symdyn(orbit)
 
 cnt = 0;
 [d,n] = size(orbit);
 s = [];
 for i = 1:n-1
   if((orbit(3,i) <= 27) & (orbit(3,i+1)> 27))
       cnt = cnt+1;
       if(orbit(2,i) < 0)
         s = [s,'A'];
       else
         s = [s,'B'];
        end
     end
 end
 if((orbit(3,n) <= 27) & (orbit(3,1)> 27))
    cnt = cnt+1;
    if(orbit(2,n) < 0)
      s = [s,'A'];
    else
      s = [s,'B'];
    end
 end