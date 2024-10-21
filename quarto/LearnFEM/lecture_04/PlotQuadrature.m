function PlotQuadrature( qp, weights, func_vals )
    num_qp = length( qp );
    verts = [];
    faces = zeros( num_qp, 5 );
    colors = zeros( num_qp, 1, 3 );
    cmap = lines(2);
    for ii = 1 : length( qp )
        verts = [ verts; 
                  [ qp(ii) + weights(ii) / 2, 0 ];
                  [ qp(ii) + weights(ii) / 2, func_vals(ii) ];
                  [ qp(ii) - weights(ii) / 2, func_vals(ii) ];
                  [ qp(ii) - weights(ii) / 2, 0 ];
                ];
        faces(ii,:) = (4*(ii-1)) + [1:4 1];
        colors(ii,1,:) = cmap( ( func_vals(ii) < 0 ) + 1, : );
    end

    ax = gca;
    old_next_plot = ax.NextPlot;
    ax.NextPlot = "add";

    P = patch( NaN, NaN, NaN );
    P.Vertices = verts;
    P.Faces = faces;
    P.CData = colors;
    P.FaceAlpha = 0.5;
    P.FaceColor = "flat";
    P.EdgeColor = "k";
    
    
    S = scatter( NaN, NaN );
    S.XData = qp;
    S.YData = func_vals;
    S.SizeData = 40;
    S.MarkerFaceColor = "k";
    S.MarkerEdgeColor = "w";

    ax.NextPlot = old_next_plot;
end