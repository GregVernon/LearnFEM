function [u, M, F, basis, d] = LaplaceEquation( basis_family, degree, domain, f, g, h, boundary_tol, penalty )
variate = symvar( f );
basis = @(deriv) diff( PolynomialBasisFunction( basis_family, degree, variate, domain ), variate, deriv );
[M, F] = ApplyGoverningEquation( basis, domain, f );
[M, F] = ApplyBoundaryConditions( M, F, basis, domain, g, h, boundary_tol, penalty );
d = M \ F;
u = symfun( transpose( d ) * basis(0), variate );
end

function [M, F] = ApplyGoverningEquation( basis, domain, distributed_load )
    M = int( basis(2) .* transpose( basis(2) ), domain );
    F = -int( basis(2) * distributed_load, domain );
end

function [M, F] = ApplyBoundaryConditions( M, F, basis, domain, g, h, boundary_tol, penalty )
[MD, FD] = DirichletBoundaryCondition( basis, domain, g, boundary_tol, penalty );
[MN, FN] = NeumannBoundaryCondition( basis, domain, h, boundary_tol, penalty );
M = M + MD + MN;
F = F + FD + FN;
end

function [M, F] = DirichletBoundaryCondition( basis, domain, g, boundary_tol, penalty )
M = penalty * int( basis(0) .* transpose( basis(0) ), [domain(2) - boundary_tol, domain(2) + boundary_tol] );
F = penalty * int( basis(0) .* g,              [domain(2) - boundary_tol, domain(2) + boundary_tol] );
end

function [M, F] = NeumannBoundaryCondition( basis, domain, h, boundary_tol, penalty )
M = penalty * int( basis(1) .* transpose( basis(1) ), [domain(1) - boundary_tol, domain(1) + boundary_tol] );
F = penalty * int( basis(1) .* h, [domain(1) - boundary_tol, domain(1) + boundary_tol] );
end