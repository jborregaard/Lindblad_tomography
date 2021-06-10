function [Out] = evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p00)
%function that calculates evolution due to Kraus operators

Out=K1*p00*ctranspose(K1)+K2*p00*ctranspose(K2)+K3*p00*ctranspose(K3)+K4*p00*ctranspose(K4)+K5*p00*ctranspose(K5)+K6*p00*ctranspose(K6)+K7*p00*ctranspose(K7)+K8*p00*ctranspose(K8)+K9*p00*ctranspose(K9)+K10*p00*ctranspose(K10)+K11*p00*ctranspose(K11)+K12*p00*ctranspose(K12)+K13*p00*ctranspose(K13)+K14*p00*ctranspose(K14)+K15*p00*ctranspose(K15)+K16*p00*ctranspose(K16);

end

