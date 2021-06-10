function [c,Out] = constraintpovm(xm)
%This function is to force the contraint that the sum of POVM elements
%gives the identity. 
c=[];

POVMo(1,1)=xm(1,1);
POVMo(2,1)=xm(1,2)+1i*xm(1,3);
POVMo(2,2)=xm(1,4);

M0=POVMo*ctranspose(POVMo);


POVM=M0;

O=eye(2)-POVM;

b=eig(O);

Out=[b(1,1)-abs(b(1,1)) b(2,1)-abs(b(2,1))];
end

