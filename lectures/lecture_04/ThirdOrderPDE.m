function [u, M, F, basis, d] = ThirdOrderPDE( basis_family, degree, domain, distributed_load, displacement_constraint, slope_constraint, gradient_constraint )
variate = symvar( distributed_load );
basis = @(deriv) diff( PolynomialBasisFunction( basis_family, degree, variate, domain ), variate, deriv );
[M, F] = ApplyGoverningEquation( basis, domain, distributed_load );
[M, F] = ApplyBoundaryConditions( M, F, basis, domain, displacement_constraint, slope_constraint, gradient_constraint );
d = M \ F;
u = symfun( transpose( d ) * basis(0), variate );
u = simplify( u, Steps=10 );

    function [M,F] = ApplyGoverningEquation( basis, domain, distributed_load )
        M = int( basis(3) .* transpose( basis(3) ), domain );
        F = -1 * int( basis(3) * distributed_load, domain );
    end

    function [M, F] = ApplyBoundaryConditions( M, F, basis, domain, displacement_constraint, slope_constraint, gradient_constraint )
        [MD, FD] = DisplacementBoundaryCondition( basis, domain, displacement_constraint );
        [MN, FN] = SlopeBoundaryCondition( basis, domain, slope_constraint );
        [MG, FG] = GradientBoundaryCondition(basis, domain, gradient_constraint );
        M = M + MD + MN + MG;
        F = F + FD + FN + FG;
    end

    function [M, F] = DisplacementBoundaryCondition( basis, domain, displacement_constraint )
        M = subs( basis(0) .* transpose( basis(0) ), symvar( basis(0) ), domain(1) );
        F = subs( basis(0) .* displacement_constraint, symvar( basis(0) ), domain(1) );
    end

    function [M, F] = SlopeBoundaryCondition( basis, domain, slope_constraint )
        M = subs( basis(1) .* transpose( basis(1) ), symvar( basis(0) ), domain(1) );
        F = subs( basis(1) .* slope_constraint, symvar( basis(0) ), domain(1) );
    end

    function [M, F] = GradientBoundaryCondition( basis, domain, gradient_constraint )
        M = subs( basis(2) .* transpose( basis(2) ), symvar( basis(0) ), domain(1) );
        F = subs( basis(2) .* gradient_constraint, symvar( basis(0) ), domain(1) );
    end
end