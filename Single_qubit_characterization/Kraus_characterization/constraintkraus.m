function [c,Out] = constraintkraus(xm)
%This function is to force the contraint that the sum of Kraus operators
%gives the identity. 
c=[];

k11=xm(1,1)+1i*xm(1,2);
k12=xm(1,3)+1i*xm(1,4);
k13=xm(1,5)+1i*xm(1,6);
k14=xm(1,7)+1i*xm(1,8);

k21=xm(1,9)+1i*xm(1,10); 
k22=xm(1,11)+1i*xm(1,12);
k23=xm(1,13)+1i*xm(1,14);
k24=xm(1,15)+1i*xm(1,16);

k31=xm(1,17)+1i*xm(1,18);
k32=xm(1,19)+1i*xm(1,20);
k33=xm(1,21)+1i*xm(1,22);
k34=xm(1,23)+1i*xm(1,24);

k41=xm(1,25)+1i*xm(1,26);
k42=xm(1,27)+1i*xm(1,28);
k43=xm(1,29)+1i*xm(1,30);
k44=xm(1,31)+1i*xm(1,32);


O(1,1)=k11*conj(k11)+k13*conj(k13)+k21*conj(k21)+k23*conj(k23)+k31*conj(k31)+k33*conj(k33)+k41*conj(k41)+k43*conj(k43);
O(1,2)=k12*conj(k11)+k14*conj(k13)+k22*conj(k21)+k24*conj(k23)+k32*conj(k31)+k34*conj(k33)+k42*conj(k41)+k44*conj(k43);
O(2,1)=k11*conj(k12)+k13*conj(k14)+k21*conj(k22)+k23*conj(k24)+k31*conj(k32)+k33*conj(k34)+k41*conj(k42)+k43*conj(k44);
O(2,2)=k12*conj(k12)+k14*conj(k14)+k22*conj(k22)+k24*conj(k24)+k32*conj(k32)+k34*conj(k34)+k42*conj(k42)+k44*conj(k44);

Out=[abs(O(1,1)-1) abs(O(2,2)-1) abs(O(1,2)) abs(O(2,1))];
end

