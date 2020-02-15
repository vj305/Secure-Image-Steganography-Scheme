function [proImageOut] = imageProcess(ImgInp,key,key1)
[n m k] = size(ImgInp);
for ind = 1 : m    
    Fkey(:,ind) = key((1+(ind-1)*n) : (((ind-1)*n)+n));
end
len = n;
bre = m;
for ind = 1 : k
    Img = ImgInp(:,:,ind);
for ind1 = 1 : len
    for ind2 = 1 : bre        
        proImage(ind1,ind2) = bitxor(Img(ind1,ind2),Fkey(ind1,ind2));        
    end
end
proImageOut(:,:,ind) = proImage(:,:,1);
end
proImageOut(1,end,1)=key1;
return;