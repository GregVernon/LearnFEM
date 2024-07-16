function [u, M, F, basis, d] = ScalarProjection( basis_family, target_fun, degree, domain )
variate = symvar( target_fun );
basis = PolynomialBasisFunction( basis_family, degree, variate, domain );
[M, F] = ApplyGoverningEquation( basis, target_fun, domain );
[M, F] = ApplyBoundaryConditions( M, F, basis, target_fun, domain );
d = M \ F;
u = symfun( transpose( d ) * basis, variate );
u = simplify( u, Steps=10 );

    function [M, F] = ApplyGoverningEquation( basis, target_fun, domain )
        M = int( basis .* transpose( basis ), domain );
        F = int( basis * target_fun, domain );
    end
    
    function [M, F] = ApplyBoundaryConditions( M, F, basis, target_fun, domain )
        [ML, FL] = LeftBoundaryCondition( basis, target_fun, domain );
        [MR, FR] = RightBoundaryCondition( basis, target_fun, domain );
        M = M + ML + MR;
        F = F + FL + FR;
    end

    function [M, F] = LeftBoundaryCondition( basis, target_fun, domain )
        M = subs( basis .* transpose( basis ), domain(1) );
        F = subs( basis * boundary_val, symvar( target_fun ), domain(1) );
    end

    function [M, F] = RightBoundaryCondition( basis, target_fun, domain )
        M = subs( basis .* transpose( basis ), domain(2) );
        F = subs( basis * boundary_val, symvar( target_fun ), domain(2) );
    end
end