classdef farea < handle

    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetObservable, AbortSet)
        Function
        XRange
        XData = NaN
        YData = NaN
        MeshDensity = 100
        EdgeColor = "k"
        FaceColor
        FaceAlpha = 0.5
        LineWidth = 0.5
        LineStyle = "-"
    end

    properties(Access = private,Transient,NonCopyable)
        AreaObject
    end

    methods( Access = public )
        function obj = farea( Function, XRange, varargin )
            p = inputParser;
            addRequired( p, "Function" );
            addRequired( p, "XRange" );
            addParameter( p, "ax", [], @(x)isgraphics(x, "matlab.graphics.axis.Axes") );
            addParameter( p, "MeshDensity", 100 );
            addParameter( p, "FaceColor", [] );
            addParameter( p, "FaceAlpha", [] );
            addParameter( p, "EdgeColor", [] );
            addParameter( p, "LineStyle", [] );
            addParameter( p, "LineWidth", [] );
            parse( p, Function, XRange, varargin{:} );

            obj.Function = Function;
            obj.XRange = XRange;

            obj = init( obj, p );
            addlistener( obj, 'Function',    'PostSet', @obj.updateFunction );
            addlistener( obj, 'XRange',      'PostSet', @obj.updateXRange );
            addlistener( obj, 'MeshDensity', 'PostSet', @obj.updateMeshDensity );
            addlistener( obj, 'EdgeColor',   'PostSet', @obj.updateEdgeColor );
            addlistener( obj, 'FaceColor',   'PostSet', @obj.updateFaceColor );
            addlistener( obj, 'LineWidth',   'PostSet', @obj.updateLineWidth );
        end
    end

    methods( Access = protected )
        function obj = init( obj, inputs )
            % Preprocess the function
            if isempty( symvar( obj.Function ) )
                obj.Function = symfun( obj.Function, sym( "x" ) );
            else
                obj.Function = symfun( obj.Function, symvar( obj.Function ) );
            end
            if inputs.Results.MeshDensity
                obj.XData = linspace( obj.XRange(1), obj.XRange(2), inputs.Results.MeshDensity );
            else
                obj.XData = linspace( obj.XRange(1), obj.XRange(2), obj.MeshDensity );
            end
            obj.YData = obj.Function( obj.XData );

            % Initialize area plot
            if ~isempty( inputs.Results.ax )
                obj.AreaObject = area( inputs.Results.ax, NaN, NaN );
            else
                obj.AreaObject = area( NaN, NaN );
            end

            % Update area XData and YData
            obj.AreaObject.XData = obj.XData;
            obj.AreaObject.YData = obj.Function( obj.XData );

            % Update face
            %%% FaceColor %%%
            if ~isempty( inputs.Results.FaceColor )
                obj.FaceColor = inputs.Results.FaceColor;
            elseif isempty( obj.EdgeColor )
                obj.FaceColor = obj.EdgeColor;
            else
                ax = gca;
                obj.FaceColor = ax.ColorOrder( ax.ColorOrderIndex, : );
            end
            obj.AreaObject.FaceColor = obj.FaceColor;
            
            %%% FaceAlpha %%%
            if ~isempty( inputs.Results.FaceAlpha )
                obj.FaceAlpha = inputs.Results.FaceAlpha;
            end
            obj.AreaObject.FaceAlpha = obj.FaceAlpha;

            % Update edge
            %%% LineStyle %%%
            if ~isempty( inputs.Results.LineStyle )
                obj.LineStyle = inputs.Results.LineStyle;
            end
            obj.AreaObject.LineStyle = obj.LineStyle;

            %%% LineWidth %%%
            if ~isempty( inputs.Results.LineWidth )
                obj.LineWidth = inputs.Results.LineWidth;
            end
            obj.AreaObject.LineWidth = obj.LineWidth;

            %%% EdgeColor %%%
            if ~isempty( inputs.Results.EdgeColor )
                obj.EdgeColor = inputs.Results.EdgeColor;
            else
                ax = gca;
                obj.EdgeColor = ax.ColorOrder( ax.ColorOrderIndex, : );
            end
            obj.AreaObject.EdgeColor = obj.EdgeColor;
        end

        function updateFunction( obj, src, event )
            % Preprocess the function
             if isempty( symvar( obj.Function ) )
                obj.Function = symfun( obj.Function, sym( "x" ) );
            else
                obj.Function = symfun( obj.Function, symvar( obj.Function ) );
             end
            obj.XData = linspace( obj.XRange(1), obj.XRange(2), obj.MeshDensity );
            obj.YData = obj.Function( obj.XData );

            % Update area XData and YData
            obj.AreaObject.XData = obj.XData;
            obj.AreaObject.YData = obj.Function( obj.XData );
        end

        function updateXRange( obj, src, event )
            % Preprocess the function
            obj.Function = symfun( obj.Function, symvar( obj.Function ) );
            obj.XData = linspace( obj.XRange(1), obj.XRange(2), obj.MeshDensity );
            obj.YData = obj.Function( obj.XData );

            % Update area XData and YData
            obj.AreaObject.XData = obj.XData;
            obj.AreaObject.YData = obj.Function( obj.XData );
        end

        function updateMeshDensity( obj, src, event )
            % Preprocess the function
            obj.Function = symfun( obj.Function, symvar( obj.Function ) );
            obj.XData = linspace( obj.XRange(1), obj.XRange(2), obj.MeshDensity );
            obj.YData = obj.Function( obj.XData );

            % Update area XData and YData
            obj.AreaObject.XData = obj.XData;
            obj.AreaObject.YData = obj.Function( obj.XData );
        end

        function updateEdgeColor( obj, src, event )
            obj.AreaObject.EdgeColor = obj.EdgeColor;
        end

        function updateFaceColor( obj, src, event )
            obj.AreaObject.FaceColor = obj.FaceColor;
        end

        function updateLineWidth( obj, src, event )
            obj.AreaObject.LineWidth = obj.LineWidth;
        end
    end
    %         function obj = untitled( fun, xinterval, varargin )
    %             %UNTITLED Construct an instance of this class
    %             %   Detailed explanation goes here
    %             p = inputParser;
    %             addRequired( p, "fun" );
    %             addRequired( p, "xinterval" );
    %             addParameter( p, "ax", [], @(x)isgraphics(x, "matlab.graphics.axis.Axes") );
    %             addParameter( p, "MeshDensity", 100 );
    %             addParameter( p, "FaceColor", [] );
    %             addParameter( p, "EdgeColor", [] );
    %             addParameter( p, "LineStyle", [] );
    %             addParameter( p, "LineWidth", [] );
    %             parse( p, fun, xinterval, varargin{:} );
    %
    %             variate = symvar( fun );
    %
    %             obj.XRange = xinterval;
    %             obj.XData = linspace( xinterval(1), xinterval(2), p.Results.MeshDensity );
    %             fun = matlabFunction( fun, Vars = variate );
    %             obj.YData = fun( obj.XData );
    %         end
end