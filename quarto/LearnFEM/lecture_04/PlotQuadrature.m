function PlotQuadrature( qp, weights, func_vals, color_style )
    num_qp = length( qp );
    verts = [];
    faces = zeros( num_qp, 5 );
    if nargin == 3 || ( isstring( color_style ) && strcomp( color_style, "default" ) )
        colors = zeros( num_qp, 1, 3 );
        cmap = lines(2);
    end
    for ii = 1 : length( qp )
        verts = [ verts; 
                  [ qp(ii) + weights(ii) / 2, 0 ];
                  [ qp(ii) + weights(ii) / 2, func_vals(ii) ];
                  [ qp(ii) - weights(ii) / 2, func_vals(ii) ];
                  [ qp(ii) - weights(ii) / 2, 0 ];
                ];
        faces(ii,:) = (4*(ii-1)) + [1:4 1];
        if nargin == 3 || ( isstring( color_style ) && strcomp( color_style, "default" ) )
            colors(ii,1,:) = cmap( ( func_vals(ii) < 0 ) + 1, : );
        end
    end

    ax = gca;
    old_next_plot = ax.NextPlot;
    ax.NextPlot = "add";

    P = patch( NaN, NaN, NaN );
    P.Vertices = verts;
    P.Faces = faces;
    if nargin == 3 || ( isstring( color_style ) && strcomp( color_style, "default" ) )
        P.CData = colors;
        P.FaceColor = "flat";
    else
        P.FaceColor = color_style;
    end
    P.FaceAlpha = 0.5;
    P.EdgeColor = "k";
    
    
    S = scatter( NaN, NaN );
    S.XData = qp;
    S.YData = func_vals;
    S.SizeData = 40;
    S.MarkerFaceColor = "k";
    S.MarkerEdgeColor = "w";

    ax.NextPlot = old_next_plot;
end