function [s,s1]=psnr(A,B)
[n,m,ch]=size(B);
A=double(A);
B=double(B);
if ch==1
   e=A-B;
   e=e(0+1:n-0,0+1:m-0);
   me=mean(mean(e.^2));
   s=10*log10(255^2/me);
   s1=me;
else
   e=A-B;
   e=e(0+1:n-0,0+1:m-0,:);
   e1=e(:,:,1);e2=e(:,:,2);e3=e(:,:,3);
   me1=mean(mean(e1.^2));
   me2=mean(mean(e2.^2));
   me3=mean(mean(e3.^2));
   mse=(me1+me2+me3)/3;
   s  = 10*log10(255^2/mse);
   s1=mse;
end


