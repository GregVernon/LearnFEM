function [u, M, F, basis, d] = PoissonEquation( spline_space, f, g, h )
variate = symvar(f);
SplineBasis = PiecewisePolynomialBasisFunction( spline_space, variate );
basis =  @(deriv) diff( SplineBasis.Basis, variate, deriv );
domain = spline_space.getSplineDomain();
[M, F] = ApplyGoverningEquation( basis, domain, f );
[M, F] = ApplyBoundaryConditions( M, F, basis, domain, g, h );
d = M \ F;
u = symfun( transpose( d ) * basis(0), variate );
u = simplify( u, Steps=10 );

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
        M = subs( basis(0) .* transpose( basis(0) ), variate, domain(2) );
        F = subs( basis(0) * g, variate, domain(2) );
    end

    function [M, F] = NeumannBoundaryCondition( basis, domain, h )
        M = limit( basis(1), variate, domain(1), "right" ) .* transpose( limit( basis(1), variate, domain(1), "right" ) );
        F = limit( basis(1), variate, domain(1), "right" ) * h;
    end
end