%Initial point for SPAM characterization

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



    
    
    
    
%%Sample rho
Q = orth(randn(2) + randn(2)*1i);
D = diag(abs(randn(2, 1)));
A = Q*D*Q';
rho0n=nearestSPD(A);
rho0=1/trace(rho0n)*rho0n;