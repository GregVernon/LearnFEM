function [u, M, F, basis, d] = ScalarProjection( spline_space, target_fun )
variate = symvar( target_fun );
SplineBasis = PiecewisePolynomialBasisFunction( spline_space, variate );
basis = SplineBasis.Basis;
domain = spline_space.getSplineDomain();
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
        F = subs( basis * target_fun, symvar( target_fun ), domain(1) );
    end

    function [M, F] = RightBoundaryCondition( basis, target_fun, domain )
        M = subs( basis .* transpose( basis ), domain(2) );
        F = subs( basis * target_fun, symvar( target_fun ), domain(2) );
    end
end