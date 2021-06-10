function [c,Out] = constraint_kraus1(Xm)
%This function is to force the contraint that the sum of Kraus operators
%gives the identity. 
c=[];

K1=reshape(Xm(1,1:16),4,4)+1i*reshape(Xm(1,256+1:256+16),4,4);
K2=reshape(Xm(1,1+16:2*16),4,4)+1i*reshape(Xm(1,256+1+16:256+2*16),4,4);
K3=reshape(Xm(1,1+2*16:3*16),4,4)+1i*reshape(Xm(1,256+1+2*16:256+3*16),4,4);
K4=reshape(Xm(1,1+3*16:4*16),4,4)+1i*reshape(Xm(1,256+1+3*16:256+4*16),4,4);
K5=reshape(Xm(1,1+4*16:5*16),4,4)+1i*reshape(Xm(1,256+1+4*16:256+5*16),4,4);
K6=reshape(Xm(1,1+5*16:6*16),4,4)+1i*reshape(Xm(1,256+1+5*16:256+6*16),4,4);
K7=reshape(Xm(1,1+6*16:7*16),4,4)+1i*reshape(Xm(1,256+1+6*16:256+7*16),4,4);
K8=reshape(Xm(1,1+7*16:8*16),4,4)+1i*reshape(Xm(1,256+1+7*16:256+8*16),4,4);
K9=reshape(Xm(1,1+8*16:9*16),4,4)+1i*reshape(Xm(1,256+1+8*16:256+9*16),4,4);
K10=reshape(Xm(1,1+9*16:10*16),4,4)+1i*reshape(Xm(1,256+1+9*16:256+10*16),4,4);
K11=reshape(Xm(1,1+10*16:11*16),4,4)+1i*reshape(Xm(1,256+1+10*16:256+11*16),4,4);
K12=reshape(Xm(1,1+11*16:12*16),4,4)+1i*reshape(Xm(1,256+1+11*16:256+12*16),4,4);
K13=reshape(Xm(1,1+12*16:13*16),4,4)+1i*reshape(Xm(1,256+1+12*16:256+13*16),4,4);
K14=reshape(Xm(1,1+13*16:14*16),4,4)+1i*reshape(Xm(1,256+1+13*16:256+14*16),4,4);
K15=reshape(Xm(1,1+14*16:15*16),4,4)+1i*reshape(Xm(1,256+1+14*16:256+15*16),4,4);
K16=reshape(Xm(1,1+15*16:16*16),4,4)+1i*reshape(Xm(1,256+1+15*16:256+16*16),4,4);

O=ctranspose(K1)*K1+ctranspose(K2)*K2+ctranspose(K3)*K3+ctranspose(K4)*K4+ctranspose(K5)*K5+ctranspose(K6)*K6+ctranspose(K7)*K7+ctranspose(K8)*K8+ctranspose(K9)*K9+ctranspose(K10)*K10+ctranspose(K11)*K11+ctranspose(K12)*K12+ctranspose(K13)*K13+ctranspose(K14)*K14+ctranspose(K15)*K15+ctranspose(K16)*K16;

Out=reshape(abs(O-eye(4,4)),1,16);
end

