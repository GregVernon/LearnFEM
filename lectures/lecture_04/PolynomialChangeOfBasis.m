function [d, R, D, C] = PolynomialChangeOfBasis( fromBasis, toBasis, fromCoeff, toDomain )
N = numel( toBasis );
D = sym( zeros( N ) );
for ii = 1 : N
    for jj = 1 : N
        D(ii,jj) = int( toBasis(ii) * toBasis(jj), toDomain );
    end
end

N1 = numel( toBasis );
N2 = numel( fromBasis );
C = sym( zeros( N1, N2 ) );
for ii = 1 : N1
    for jj = 1 : N2
        C(ii,jj) = int( toBasis(ii) * fromBasis(jj), toDomain );
    end
end

R = D \ C;
d = R * fromCoeff;
end