function poly_roots = PolynomialRoots( polyfun )
C = CompanionMatrix( polyfun );
poly_roots = eig( C );
[~, sidx] = sort( double ( real( vpa( poly_roots ) ) ) );
poly_roots = simplify( poly_roots( sidx ), Steps=100, Seconds=1 ); 
end