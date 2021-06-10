function [Out] = likelihood_povm2(f,M,rho,c)
%Function with the SPAM likelihood functions

%The # of measurements for all 36 input states. I.e. f11 is the number of
%measurements for inpout state |0>|0>, f12 is for |0>|1>, f13 is for
%|0>|+>, f14 is for |0>|->, f15 is for |0>|+i>, f16 is for |0>|-i>, f21 is
%for |1>|0> etc. Note that the measurment basis and outcome has already
%been specified in likelihood_povm1. 
f11=f(1,1); f12=f(1,2); f13=f(1,3); f14=f(1,4); f15=f(1,5); f16=f(1,6);
f21=f(2,1); f22=f(2,2); f23=f(2,3); f24=f(2,4); f25=f(2,5); f26=f(2,6);
f31=f(3,1); f32=f(3,2); f33=f(3,3); f34=f(3,4); f35=f(3,5); f36=f(3,6);
f41=f(4,1); f42=f(4,2); f43=f(4,3); f44=f(4,4); f45=f(4,5); f46=f(4,6);
f51=f(5,1); f52=f(5,2); f53=f(5,3); f54=f(5,4); f55=f(5,5); f56=f(5,6);
f61=f(6,1); f62=f(6,2); f63=f(6,3); f64=f(6,4); f65=f(6,5); f66=f(6,6);

%the single qubit Gates for initial state preparation:
sigmax=[0 -1i;-1i 0];
hadp=1/sqrt(2)*[1 -1;1 1];
hadip=1/sqrt(2)*[1 1i; 1i 1];
hadm=1/sqrt(2)*[1 -1;1 1];
hadim=1/sqrt(2)*[1 -1i; -1i 1];



%The 36 different initial two qubit states
p00=rho;
p01=kron(eye(2),sigmax)*rho*ctranspose(kron(eye(2),sigmax));
p0p=kron(eye(2),hadp)*rho*ctranspose(kron(eye(2),hadp));
p0m=kron(eye(2),hadm)*rho*ctranspose(kron(eye(2),hadm));
p0ip=kron(eye(2),hadip)*rho*ctranspose(kron(eye(2),hadip));
p0im=kron(eye(2),hadim)*rho*ctranspose(kron(eye(2),hadim));

p10=kron(sigmax,eye(2))*rho*ctranspose(kron(sigmax,eye(2)));
p11=kron(sigmax,sigmax)*rho*ctranspose(kron(sigmax,sigmax));
p1p=kron(sigmax,hadp)*rho*ctranspose(kron(sigmax,hadp));
p1m=kron(sigmax,hadm)*rho*ctranspose(kron(sigmax,hadm));
p1ip=kron(sigmax,hadip)*rho*ctranspose(kron(sigmax,hadip));
p1im=kron(sigmax,hadim)*rho*ctranspose(kron(sigmax,hadim));

pp0=kron(hadp,eye(2))*rho*ctranspose(kron(hadp,eye(2)));
pp1=kron(hadp,sigmax)*rho*ctranspose(kron(hadp,sigmax));
ppp=kron(hadp,hadp)*rho*ctranspose(kron(hadp,hadp));
ppm=kron(hadp,hadm)*rho*ctranspose(kron(hadp,hadm));
ppip=kron(hadp,hadip)*rho*ctranspose(kron(hadp,hadip));
ppim=kron(hadp,hadim)*rho*ctranspose(kron(hadp,hadim));

pm0=kron(hadm,eye(2))*rho*ctranspose(kron(hadm,eye(2)));
pm1=kron(hadm,sigmax)*rho*ctranspose(kron(hadm,sigmax));
pmp=kron(hadm,hadp)*rho*ctranspose(kron(hadm,hadp));
pmm=kron(hadm,hadm)*rho*ctranspose(kron(hadm,hadm));
pmip=kron(hadm,hadip)*rho*ctranspose(kron(hadm,hadip));
pmim=kron(hadm,hadim)*rho*ctranspose(kron(hadm,hadim));

pip0=kron(hadip,eye(2))*rho*ctranspose(kron(hadip,eye(2)));
pip1=kron(hadip,sigmax)*rho*ctranspose(kron(hadip,sigmax));
pipp=kron(hadip,hadp)*rho*ctranspose(kron(hadip,hadp));
pipm=kron(hadip,hadm)*rho*ctranspose(kron(hadip,hadm));
pipip=kron(hadip,hadip)*rho*ctranspose(kron(hadip,hadip));
pipim=kron(hadip,hadim)*rho*ctranspose(kron(hadip,hadim));

pim0=kron(hadim,eye(2))*rho*ctranspose(kron(hadim,eye(2)));
pim1=kron(hadim,sigmax)*rho*ctranspose(kron(hadim,sigmax));
pimp=kron(hadim,hadp)*rho*ctranspose(kron(hadim,hadp));
pimm=kron(hadim,hadm)*rho*ctranspose(kron(hadim,hadm));
pimip=kron(hadim,hadip)*rho*ctranspose(kron(hadim,hadip));
pimim=kron(hadim,hadim)*rho*ctranspose(kron(hadim,hadim));

if c==1
%Lzz
Mb=M;
    Out=f11*log(trace(Mb*p00))+f12*log(trace(Mb*p01))+f13*log(trace(Mb*p0p))+f14*log(trace(Mb*p0m))+f15*log(trace(Mb*p0ip))+f16*log(trace(Mb*p0im))+...
        f21*log(trace(Mb*p10))+f22*log(trace(Mb*p11))+f23*log(trace(Mb*p1p))+f24*log(trace(Mb*p1m))+f25*log(trace(Mb*p1ip))+f26*log(trace(Mb*p1im))+...
        f31*log(trace(Mb*pp0))+f32*log(trace(Mb*pp1))+f33*log(trace(Mb*ppp))+f34*log(trace(Mb*ppm))+f35*log(trace(Mb*ppip))+f36*log(trace(Mb*ppim))+...
        f41*log(trace(Mb*pm0))+f42*log(trace(Mb*pm1))+f43*log(trace(Mb*pmp))+f44*log(trace(Mb*pmm))+f45*log(trace(Mb*pmip))+f46*log(trace(Mb*pmim))+...
        f51*log(trace(Mb*pip0))+f52*log(trace(Mb*pip1))+f53*log(trace(Mb*pipp))+f54*log(trace(Mb*pipm))+f55*log(trace(Mb*pipip))+f56*log(trace(Mb*pipim))+...
        f61*log(trace(Mb*pim0))+f62*log(trace(Mb*pim1))+f63*log(trace(Mb*pimp))+f64*log(trace(Mb*pimm))+f65*log(trace(Mb*pimip))+f66*log(trace(Mb*pimim));
end
if c==2
%Lzx
Mzx=[2.^(-1/2),2.^(-1/2),0,0;(-1).*2.^(-1/2),2.^(-1/2),0,0;0,0,2.^(-1/2),2.^(-1/2);0,0,(-1).*2.^(-1/2),2.^(-1/2)];
Mb=ctranspose(Mzx)*M*Mzx;
    Out=f11*log(trace(Mb*p00))+f12*log(trace(Mb*p01))+f13*log(trace(Mb*p0p))+f14*log(trace(Mb*p0m))+f15*log(trace(Mb*p0ip))+f16*log(trace(Mb*p0im))+...
        f21*log(trace(Mb*p10))+f22*log(trace(Mb*p11))+f23*log(trace(Mb*p1p))+f24*log(trace(Mb*p1m))+f25*log(trace(Mb*p1ip))+f26*log(trace(Mb*p1im))+...
        f31*log(trace(Mb*pp0))+f32*log(trace(Mb*pp1))+f33*log(trace(Mb*ppp))+f34*log(trace(Mb*ppm))+f35*log(trace(Mb*ppip))+f36*log(trace(Mb*ppim))+...
        f41*log(trace(Mb*pm0))+f42*log(trace(Mb*pm1))+f43*log(trace(Mb*pmp))+f44*log(trace(Mb*pmm))+f45*log(trace(Mb*pmip))+f46*log(trace(Mb*pmim))+...
        f51*log(trace(Mb*pip0))+f52*log(trace(Mb*pip1))+f53*log(trace(Mb*pipp))+f54*log(trace(Mb*pipm))+f55*log(trace(Mb*pipip))+f56*log(trace(Mb*pipim))+...
        f61*log(trace(Mb*pim0))+f62*log(trace(Mb*pim1))+f63*log(trace(Mb*pimp))+f64*log(trace(Mb*pimm))+f65*log(trace(Mb*pimip))+f66*log(trace(Mb*pimim));
end
if c==3
%Lzy
Mzy=[2.^(-1/2),(sqrt(-1)*(-1)).*2.^(-1/2),0,0;(sqrt(-1)*(-1)).*2.^(-1/2),2.^(-1/2),0,0;0,0,2.^(-1/2),(sqrt(-1)*(-1)).*2.^(-1/2);0,0,(sqrt(-1)*(-1)).*2.^(-1/2),2.^(-1/2)];
Mb=ctranspose(Mzy)*M*Mzy;
    Out=f11*log(trace(Mb*p00))+f12*log(trace(Mb*p01))+f13*log(trace(Mb*p0p))+f14*log(trace(Mb*p0m))+f15*log(trace(Mb*p0ip))+f16*log(trace(Mb*p0im))+...
        f21*log(trace(Mb*p10))+f22*log(trace(Mb*p11))+f23*log(trace(Mb*p1p))+f24*log(trace(Mb*p1m))+f25*log(trace(Mb*p1ip))+f26*log(trace(Mb*p1im))+...
        f31*log(trace(Mb*pp0))+f32*log(trace(Mb*pp1))+f33*log(trace(Mb*ppp))+f34*log(trace(Mb*ppm))+f35*log(trace(Mb*ppip))+f36*log(trace(Mb*ppim))+...
        f41*log(trace(Mb*pm0))+f42*log(trace(Mb*pm1))+f43*log(trace(Mb*pmp))+f44*log(trace(Mb*pmm))+f45*log(trace(Mb*pmip))+f46*log(trace(Mb*pmim))+...
        f51*log(trace(Mb*pip0))+f52*log(trace(Mb*pip1))+f53*log(trace(Mb*pipp))+f54*log(trace(Mb*pipm))+f55*log(trace(Mb*pipip))+f56*log(trace(Mb*pipim))+...
        f61*log(trace(Mb*pim0))+f62*log(trace(Mb*pim1))+f63*log(trace(Mb*pimp))+f64*log(trace(Mb*pimm))+f65*log(trace(Mb*pimip))+f66*log(trace(Mb*pimim));
end

if c==4
%Lxz
Mxz=[2.^(-1/2),0,2.^(-1/2),0;0,2.^(-1/2),0,2.^(-1/2);(-1).*2.^(-1/2),0,2.^(-1/2),0;0,(-1).*2.^(-1/2),0,2.^(-1/2)];
Mb=ctranspose(Mxz)*M*Mxz;
    Out=f11*log(trace(Mb*p00))+f12*log(trace(Mb*p01))+f13*log(trace(Mb*p0p))+f14*log(trace(Mb*p0m))+f15*log(trace(Mb*p0ip))+f16*log(trace(Mb*p0im))+...
        f21*log(trace(Mb*p10))+f22*log(trace(Mb*p11))+f23*log(trace(Mb*p1p))+f24*log(trace(Mb*p1m))+f25*log(trace(Mb*p1ip))+f26*log(trace(Mb*p1im))+...
        f31*log(trace(Mb*pp0))+f32*log(trace(Mb*pp1))+f33*log(trace(Mb*ppp))+f34*log(trace(Mb*ppm))+f35*log(trace(Mb*ppip))+f36*log(trace(Mb*ppim))+...
        f41*log(trace(Mb*pm0))+f42*log(trace(Mb*pm1))+f43*log(trace(Mb*pmp))+f44*log(trace(Mb*pmm))+f45*log(trace(Mb*pmip))+f46*log(trace(Mb*pmim))+...
        f51*log(trace(Mb*pip0))+f52*log(trace(Mb*pip1))+f53*log(trace(Mb*pipp))+f54*log(trace(Mb*pipm))+f55*log(trace(Mb*pipip))+f56*log(trace(Mb*pipim))+...
        f61*log(trace(Mb*pim0))+f62*log(trace(Mb*pim1))+f63*log(trace(Mb*pimp))+f64*log(trace(Mb*pimm))+f65*log(trace(Mb*pimip))+f66*log(trace(Mb*pimim));
end
if c==5
%Lxx
Mxx=[(1/2),(1/2),(1/2),(1/2);(-1/2),(1/2),(-1/2),(1/2);(-1/2),(-1/2),(1/2),(1/2);(1/2),(-1/2),(-1/2),(1/2)];
Mb=ctranspose(Mxx)*M*Mxx;
    Out=f11*log(trace(Mb*p00))+f12*log(trace(Mb*p01))+f13*log(trace(Mb*p0p))+f14*log(trace(Mb*p0m))+f15*log(trace(Mb*p0ip))+f16*log(trace(Mb*p0im))+...
        f21*log(trace(Mb*p10))+f22*log(trace(Mb*p11))+f23*log(trace(Mb*p1p))+f24*log(trace(Mb*p1m))+f25*log(trace(Mb*p1ip))+f26*log(trace(Mb*p1im))+...
        f31*log(trace(Mb*pp0))+f32*log(trace(Mb*pp1))+f33*log(trace(Mb*ppp))+f34*log(trace(Mb*ppm))+f35*log(trace(Mb*ppip))+f36*log(trace(Mb*ppim))+...
        f41*log(trace(Mb*pm0))+f42*log(trace(Mb*pm1))+f43*log(trace(Mb*pmp))+f44*log(trace(Mb*pmm))+f45*log(trace(Mb*pmip))+f46*log(trace(Mb*pmim))+...
        f51*log(trace(Mb*pip0))+f52*log(trace(Mb*pip1))+f53*log(trace(Mb*pipp))+f54*log(trace(Mb*pipm))+f55*log(trace(Mb*pipip))+f56*log(trace(Mb*pipim))+...
        f61*log(trace(Mb*pim0))+f62*log(trace(Mb*pim1))+f63*log(trace(Mb*pimp))+f64*log(trace(Mb*pimm))+f65*log(trace(Mb*pimip))+f66*log(trace(Mb*pimim));
end
if c==6
%Lxy
Mxy=[(1/2),(sqrt(-1)*(-1/2)),(1/2),(sqrt(-1)*(-1/2));(sqrt(-1)*(-1/2)),(1/2),(sqrt(-1)*(-1/2)),(1/2);(-1/2),(sqrt(-1)*(1/2)),(1/2),(sqrt(-1)*(-1/2));(sqrt(-1)*(1/2)),(-1/2),(sqrt(-1)*(-1/2)),(1/2)];
Mb=ctranspose(Mxy)*M*Mxy;
    Out=f11*log(trace(Mb*p00))+f12*log(trace(Mb*p01))+f13*log(trace(Mb*p0p))+f14*log(trace(Mb*p0m))+f15*log(trace(Mb*p0ip))+f16*log(trace(Mb*p0im))+...
        f21*log(trace(Mb*p10))+f22*log(trace(Mb*p11))+f23*log(trace(Mb*p1p))+f24*log(trace(Mb*p1m))+f25*log(trace(Mb*p1ip))+f26*log(trace(Mb*p1im))+...
        f31*log(trace(Mb*pp0))+f32*log(trace(Mb*pp1))+f33*log(trace(Mb*ppp))+f34*log(trace(Mb*ppm))+f35*log(trace(Mb*ppip))+f36*log(trace(Mb*ppim))+...
        f41*log(trace(Mb*pm0))+f42*log(trace(Mb*pm1))+f43*log(trace(Mb*pmp))+f44*log(trace(Mb*pmm))+f45*log(trace(Mb*pmip))+f46*log(trace(Mb*pmim))+...
        f51*log(trace(Mb*pip0))+f52*log(trace(Mb*pip1))+f53*log(trace(Mb*pipp))+f54*log(trace(Mb*pipm))+f55*log(trace(Mb*pipip))+f56*log(trace(Mb*pipim))+...
        f61*log(trace(Mb*pim0))+f62*log(trace(Mb*pim1))+f63*log(trace(Mb*pimp))+f64*log(trace(Mb*pimm))+f65*log(trace(Mb*pimip))+f66*log(trace(Mb*pimim));
end

if c==7
%Lyz
Myz=[2.^(-1/2),0,(sqrt(-1)*(-1)).*2.^(-1/2),0;0,2.^(-1/2),0,(sqrt(-1)*(-1)).*2.^(-1/2);(sqrt(-1)*(-1)).*2.^(-1/2),0,2.^(-1/2),0;0,(sqrt(-1)*(-1)).*2.^(-1/2),0,2.^(-1/2)];
Mb=ctranspose(Myz)*M*Myz;
    Out=f11*log(trace(Mb*p00))+f12*log(trace(Mb*p01))+f13*log(trace(Mb*p0p))+f14*log(trace(Mb*p0m))+f15*log(trace(Mb*p0ip))+f16*log(trace(Mb*p0im))+...
        f21*log(trace(Mb*p10))+f22*log(trace(Mb*p11))+f23*log(trace(Mb*p1p))+f24*log(trace(Mb*p1m))+f25*log(trace(Mb*p1ip))+f26*log(trace(Mb*p1im))+...
        f31*log(trace(Mb*pp0))+f32*log(trace(Mb*pp1))+f33*log(trace(Mb*ppp))+f34*log(trace(Mb*ppm))+f35*log(trace(Mb*ppip))+f36*log(trace(Mb*ppim))+...
        f41*log(trace(Mb*pm0))+f42*log(trace(Mb*pm1))+f43*log(trace(Mb*pmp))+f44*log(trace(Mb*pmm))+f45*log(trace(Mb*pmip))+f46*log(trace(Mb*pmim))+...
        f51*log(trace(Mb*pip0))+f52*log(trace(Mb*pip1))+f53*log(trace(Mb*pipp))+f54*log(trace(Mb*pipm))+f55*log(trace(Mb*pipip))+f56*log(trace(Mb*pipim))+...
        f61*log(trace(Mb*pim0))+f62*log(trace(Mb*pim1))+f63*log(trace(Mb*pimp))+f64*log(trace(Mb*pimm))+f65*log(trace(Mb*pimip))+f66*log(trace(Mb*pimim));
end
if c==8
%Lyx
Myx=[(1/2),(1/2),(sqrt(-1)*(-1/2)),(sqrt(-1)*(-1/2));(-1/2),(1/2),(sqrt(-1)*(1/2)),(sqrt(-1)*(-1/2));(sqrt(-1)*(-1/2)),(sqrt(-1)*(-1/2)),(1/2),(1/2);(sqrt(-1)*(1/2)),(sqrt(-1)*(-1/2)),(-1/2),(1/2)];
Mb=ctranspose(Myx)*M*Myx;
    Out=f11*log(trace(Mb*p00))+f12*log(trace(Mb*p01))+f13*log(trace(Mb*p0p))+f14*log(trace(Mb*p0m))+f15*log(trace(Mb*p0ip))+f16*log(trace(Mb*p0im))+...
        f21*log(trace(Mb*p10))+f22*log(trace(Mb*p11))+f23*log(trace(Mb*p1p))+f24*log(trace(Mb*p1m))+f25*log(trace(Mb*p1ip))+f26*log(trace(Mb*p1im))+...
        f31*log(trace(Mb*pp0))+f32*log(trace(Mb*pp1))+f33*log(trace(Mb*ppp))+f34*log(trace(Mb*ppm))+f35*log(trace(Mb*ppip))+f36*log(trace(Mb*ppim))+...
        f41*log(trace(Mb*pm0))+f42*log(trace(Mb*pm1))+f43*log(trace(Mb*pmp))+f44*log(trace(Mb*pmm))+f45*log(trace(Mb*pmip))+f46*log(trace(Mb*pmim))+...
        f51*log(trace(Mb*pip0))+f52*log(trace(Mb*pip1))+f53*log(trace(Mb*pipp))+f54*log(trace(Mb*pipm))+f55*log(trace(Mb*pipip))+f56*log(trace(Mb*pipim))+...
        f61*log(trace(Mb*pim0))+f62*log(trace(Mb*pim1))+f63*log(trace(Mb*pimp))+f64*log(trace(Mb*pimm))+f65*log(trace(Mb*pimip))+f66*log(trace(Mb*pimim));
end
if c==9
%Lyy
Myy=[(1/2),(sqrt(-1)*(-1/2)),(sqrt(-1)*(-1/2)),(-1/2);(sqrt(-1)*(-1/2)),(1/2),(-1/2),(sqrt(-1)*(-1/2));(sqrt(-1)*(-1/2)),(-1/2),(1/2),(sqrt(-1)*(-1/2));(-1/2),(sqrt(-1)*(-1/2)),(sqrt(-1)*(-1/2)),(1/2)];
Mb=ctranspose(Myy)*M*Myy;
    Out=f11*log(trace(Mb*p00))+f12*log(trace(Mb*p01))+f13*log(trace(Mb*p0p))+f14*log(trace(Mb*p0m))+f15*log(trace(Mb*p0ip))+f16*log(trace(Mb*p0im))+...
        f21*log(trace(Mb*p10))+f22*log(trace(Mb*p11))+f23*log(trace(Mb*p1p))+f24*log(trace(Mb*p1m))+f25*log(trace(Mb*p1ip))+f26*log(trace(Mb*p1im))+...
        f31*log(trace(Mb*pp0))+f32*log(trace(Mb*pp1))+f33*log(trace(Mb*ppp))+f34*log(trace(Mb*ppm))+f35*log(trace(Mb*ppip))+f36*log(trace(Mb*ppim))+...
        f41*log(trace(Mb*pm0))+f42*log(trace(Mb*pm1))+f43*log(trace(Mb*pmp))+f44*log(trace(Mb*pmm))+f45*log(trace(Mb*pmip))+f46*log(trace(Mb*pmim))+...
        f51*log(trace(Mb*pip0))+f52*log(trace(Mb*pip1))+f53*log(trace(Mb*pipp))+f54*log(trace(Mb*pipm))+f55*log(trace(Mb*pipip))+f56*log(trace(Mb*pipim))+...
        f61*log(trace(Mb*pim0))+f62*log(trace(Mb*pim1))+f63*log(trace(Mb*pimp))+f64*log(trace(Mb*pimm))+f65*log(trace(Mb*pimip))+f66*log(trace(Mb*pimim));
end

%Due to limited numerical precision matlab can give NaN or InF numbers
%which would kill the optimization. The fix here is to give a large penalty
%for such points but keep the optimization running. A warning is at the
%same time being issued. 
if isnan(Out)
Out=-sum(sum(f));
disp('NaN in log-likelihood');
end
if isinf(Out)
Out=-sum(sum(f));
disp('InF in log-likelihood');
end

end

