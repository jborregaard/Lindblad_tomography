function [Out] = sample_POVM(m0z,m1z,mpz,mipz)
%Function that samples a random single qubit POVM by adding "noise" to an
%initial guess of an POVM

%We use the measurement results for |0>, |1>, |+> and |+i> in the z-basis to get an
%initial guess of the POVM: 
p0=1-sum(m0z)/length(m0z);
p1=1-sum(m1z)/length(m1z);
pp=1-sum(mpz)/length(mpz);
pip=1-sum(mipz)/length(mipz);

POVM0p(1,1)=p0;
POVM0p(1,2)=real(pp-1/2*(p0+p1)-1i.*(pip-1/2*(p0+p1)))+1i*imag(pp-1/2*(p0+p1)-1i.*(pip-1/2*(p0+p1)));
POVM0p(2,1)=conj(POVM0p(1,2));
POVM0p(2,2)=p1;

%%Adding "noise"
M0=POVM0p;
M10=eye(2)-M00;
c=eig(M10);
if c-abs(c)
end

