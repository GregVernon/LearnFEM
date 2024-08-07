function [u, M, F, basis, d] = HeatEquation( spline_space, thermal_conductivity, q, g, h )
variate = symvar(q);
SplineBasis = PiecewisePolynomialBasisFunction( spline_space, variate );
basis =  @(deriv) diff( SplineBasis.Basis, variate, deriv );
domain = spline_space.getSplineDomain();
[M, F] = ApplyGoverningEquation( basis, domain, thermal_conductivity, q );
[M, F] = ApplyBoundaryConditions( M, F, basis, domain, g, h );
d = M \ F;
u = symfun( transpose( d ) * basis(0), variate );
u = simplify( u, Steps=10 );

    function [M, F] = ApplyGoverningEquation( basis, domain, thermal_conductivity, distributed_load )
        M = int( ( thermal_conductivity * basis(2) ) .* ( thermal_conductivity * transpose( basis(2) ) ), domain );
        F = -int( ( thermal_conductivity * basis(2) ) * distributed_load, domain );
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