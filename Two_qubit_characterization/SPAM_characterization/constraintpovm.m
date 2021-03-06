function [c,Out] = constraintpovm(Xm)
%This function is to force the contraint that the sum of POVM elemenbts
%gives the identity. 
c=[];

Xm00=Xm(1,1:16);
Xm01=Xm(1,17:32);
Xm10=Xm(1,33:48);
Xm11=Xm(1,49:64);

%The POVM's

L00(1,1)=abs(Xm00(1,1));
L00(2,2)=abs(Xm00(1,2));
L00(3,3)=abs(Xm00(1,3));
L00(4,4)=abs(Xm00(1,4));
L00(2,1)=Xm00(1,5)+1i*Xm00(1,6);
L00(3,2)=Xm00(1,7)+1i*Xm00(1,8);
L00(4,3)=Xm00(1,9)+1i*Xm00(1,10);
L00(3,1)=Xm00(1,11)+1i*Xm00(1,12);
L00(4,2)=Xm00(1,13)+1i*Xm00(1,14);
L00(4,1)=Xm00(1,15)+1i*Xm00(1,16);

M00=L00*ctranspose(L00);

L01(1,1)=abs(Xm01(1,1));
L01(2,2)=abs(Xm01(1,2));
L01(3,3)=abs(Xm01(1,3));
L01(4,4)=abs(Xm01(1,4));
L01(2,1)=Xm01(1,5)+1i*Xm01(1,6);
L01(3,2)=Xm01(1,7)+1i*Xm01(1,8);
L01(4,3)=Xm01(1,9)+1i*Xm01(1,10);
L01(3,1)=Xm01(1,11)+1i*Xm01(1,12);
L01(4,2)=Xm01(1,13)+1i*Xm01(1,14);
L01(4,1)=Xm01(1,15)+1i*Xm01(1,16);

M01=L01*ctranspose(L01);

L10(1,1)=abs(Xm10(1,1));
L10(2,2)=abs(Xm10(1,2));
L10(3,3)=abs(Xm10(1,3));
L10(4,4)=abs(Xm10(1,4));
L10(2,1)=Xm10(1,5)+1i*Xm10(1,6);
L10(3,2)=Xm10(1,7)+1i*Xm10(1,8);
L10(4,3)=Xm10(1,9)+1i*Xm10(1,10);
L10(3,1)=Xm10(1,11)+1i*Xm10(1,12);
L10(4,2)=Xm10(1,13)+1i*Xm10(1,14);
L10(4,1)=Xm10(1,15)+1i*Xm10(1,16);

M10=L10*ctranspose(L10);

L11(1,1)=abs(Xm11(1,1));
L11(2,2)=abs(Xm11(1,2));
L11(3,3)=abs(Xm11(1,3));
L11(4,4)=abs(Xm11(1,4));
L11(2,1)=Xm11(1,5)+1i*Xm11(1,6);
L11(3,2)=Xm11(1,7)+1i*Xm11(1,8);
L11(4,3)=Xm11(1,9)+1i*Xm11(1,10);
L11(3,1)=Xm11(1,11)+1i*Xm11(1,12);
L11(4,2)=Xm11(1,13)+1i*Xm11(1,14);
L11(4,1)=Xm11(1,15)+1i*Xm11(1,16);

M11=L11*ctranspose(L11);

O=M00+M01+M10+M11;

Out=[abs(O(1,1)-1) abs(O(2,2)-1) abs(O(3,3)-1) abs(O(4,4)-1) abs(O(1,2)) abs(O(1,3)) abs(O(1,4)) abs(O(2,1)) abs(O(2,3)) abs(O(2,4)) abs(O(3,1)) abs(O(3,2)) abs(O(3,4)) abs(O(4,1)) abs(O(4,2)) abs(O(4,3))];
end

