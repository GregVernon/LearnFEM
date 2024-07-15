function [u, M, F, basis, d] = LaplaceEquation( basis_family, degree, domain, f, g, h )
variate = symvar( f );
basis = @(deriv) diff( PolynomialBasisFunction( basis_family, degree, variate, domain ), variate, deriv );
[M, F] = ApplyGoverningEquation( basis, domain, f );
[M, F] = ApplyBoundaryConditions( M, F, basis, domain, g, h );
d = M \ F;
u = symfun( transpose( d ) * basis(0), variate );
end

function [M, F] = ApplyGoverningEquation( basis, domain, distributed_load )
    M = int( basis(2) .* transpose( basis(2) ), domain );
    F = -int( basis(2) * distributed_load, domain );
end

function [M, F] = ApplyBoundaryConditions( M, F, basis, domain, g, h )
[MD, FD] = DirichletBoundaryCondition( basis, domain, g );
[MN, FN] = NeumannBoundaryCondition( basis, domain, h );
M = M + MD + MN;
F = F + FD + FN;
end

function [M, F] = DirichletBoundaryCondition( basis, domain, g )
M = subs( basis(0) .* transpose( basis(0) ), symvar( basis(0) ), domain(2) );
F = subs( basis(0) * g, symvar( basis(0) ), domain(2) );
end

function [M, F] = NeumannBoundaryCondition( basis, domain, h )
M = subs( basis(1) .* transpose( basis(1) ), symvar( basis(1) ), domain(1) );
F = subs( basis(1) * h, symvar( basis(1) ), domain(1) );
end