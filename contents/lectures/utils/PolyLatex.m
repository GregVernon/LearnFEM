function latex_string = PolyLatex( polynomial, include_non_zero_terms, order, show_coeff_one )
% PolyLatex  Convert a univariate polynomial to a LaTeX string.
%
% This is intended for lecture notes where it's useful to:
%   - always show the variate (including the x^0 term),
%   - optionally keep / drop zero-coefficient terms,
%   - choose ascending/descending power order,
%   - optionally show coefficients of 1.
%
% Usage:
%   x = sym("x");
%   p = (1/2)*x^4 + (1/2)*x^3 - 6*x^2 - 2*x + 6;
%   PolyLatex( p )
%
% Arguments:
%   polynomial             Polynomial expression as sym, symfun, string,
%                          char vector, or numeric value. Symbolic arrays
%                          are emitted as LaTeX matrices.
%   include_non_zero_terms (logical, default=false)
%     If true, omit zero-coefficient terms.
%   order                  (string, default="descend") "descend" or "ascend".
%   show_coeff_one         (logical, default=true)
%     If false, omit coefficients of +/-1 (but still show the variate).
%
% Output:
%   latex_string (string)  LaTeX for the polynomial (e.g., "1x^2 + 0x^1 + 3x^0").

  if nargin < 2 || isempty( include_non_zero_terms )
    include_non_zero_terms = false;
  end

  if nargin < 3 || isempty( order )
    order = "descend";
  end

  if nargin < 4 || isempty( show_coeff_one )
    show_coeff_one = true;
  end

  polynomial = PolyLatexNormalizePolynomial_( polynomial );
  include_non_zero_terms = PolyLatexNormalizeLogical_( include_non_zero_terms, "include_non_zero_terms" );
  order = PolyLatexNormalizeOrder_( order );
  show_coeff_one = PolyLatexNormalizeLogical_( show_coeff_one, "show_coeff_one" );

  variate = PolyLatexGetVariate_( polynomial );
  if isscalar( polynomial )
    latex_string = PolyLatexFormatScalar_( polynomial, variate, include_non_zero_terms, order, show_coeff_one );
  else
    latex_string = PolyLatexFormatArray_( polynomial, variate, include_non_zero_terms, order, show_coeff_one );
  end

end

function variate = PolyLatexGetVariate_( polynomial )
  variate = symvar( polynomial );
  if isempty( variate )
    variate = sym( "x" );
    return
  end
  if length( variate ) > 1
    error( "PolyLatex only supports univariate polynomials." );
  end
  variate = variate(1);
end

function polynomial = PolyLatexNormalizePolynomial_( polynomial )
  if isa( polynomial, "symfun" )
    polynomial = formula( polynomial );
  elseif ischar( polynomial ) || ( isstring( polynomial ) && isscalar( polynomial ) )
    polynomial = str2sym( polynomial );
  elseif isnumeric( polynomial ) || islogical( polynomial )
    if ~isscalar( polynomial )
      error( "PolyLatex polynomial must be scalar." );
    end
    polynomial = sym( polynomial );
  end

  if ~isa( polynomial, "sym" )
    error( "PolyLatex polynomial must be sym, symfun, string, char vector, or numeric." );
  end
end

function value = PolyLatexNormalizeLogical_( value, argument_name )
  if isa( value, "logical" ) && isscalar( value )
    return
  end

  if isnumeric( value ) && isscalar( value ) && any( value == [0, 1] )
    value = logical( value );
    return
  end

  error( "PolyLatex %s must be a logical scalar (or numeric 0/1).", argument_name );
end

function order = PolyLatexNormalizeOrder_( order )
  if ~( ischar( order ) || ( isstring( order ) && isscalar( order ) ) )
    error( "PolyLatex order must be ""descend"" or ""ascend""." );
  end

  order = string( order );
  if ~ismember( order, ["descend", "ascend"] )
    error( "PolyLatex order must be ""descend"" or ""ascend""." );
  end
end

function latex_string = PolyLatexFormatScalar_( polynomial, variate, include_non_zero_terms, order, show_coeff_one )
  variate_latex = string( latex( variate ) );

  polynomial = expand( polynomial );
  c_desc = coeffs( polynomial, variate, "all" );
  c_desc = reshape( c_desc, 1, [] );
  degree = length( c_desc ) - 1;
  pow_desc = degree : -1 : 0;

  if order == "ascend"
    c = fliplr( c_desc );
    pow = fliplr( pow_desc );
  else
    c = c_desc;
    pow = pow_desc;
  end

  term = strings( 0, 1 );
  is_negative = false( 0, 1 );

  for ii = 1 : length( c )
    coeff = c(ii);
    power = pow(ii);

    if include_non_zero_terms && PolyLatexIsZero_( coeff )
      continue
    end

    [negative_coeff, abs_coeff_latex] = PolyLatexSignLatex_( coeff );

    omit_coeff = ( ~show_coeff_one ) && PolyLatexIsOneAbs_( coeff );
    if omit_coeff
      core = variate_latex + "^" + string( power );
    else
      core = abs_coeff_latex + variate_latex + "^" + string( power );
    end

    term(end+1,1) = core; %#ok<AGROW>
    is_negative(end+1,1) = negative_coeff; %#ok<AGROW>
  end

  if isempty( term )
    latex_string = "0" + variate_latex + "^0";
    return
  end

  latex_string = term(1);
  if is_negative(1)
    latex_string = "-" + latex_string;
  end

  for ii = 2 : length( term )
    if is_negative(ii)
      latex_string = latex_string + " - " + term(ii);
    else
      latex_string = latex_string + " + " + term(ii);
    end
  end
end

function latex_string = PolyLatexFormatArray_( polynomial, variate, include_non_zero_terms, order, show_coeff_one )
  [num_rows, num_cols] = size( polynomial );
  row_strings = strings( num_rows, 1 );

  for row = 1 : num_rows
    cell_strings = strings( 1, num_cols );
    for col = 1 : num_cols
      cell_strings(col) = PolyLatexFormatScalar_( polynomial(row, col), variate, include_non_zero_terms, order, show_coeff_one );
    end
    row_strings(row) = strjoin( cell_strings, " & " );
  end

  latex_string = "\begin{bmatrix}" + strjoin( row_strings, " \\ " ) + "\end{bmatrix}";
end

function tf = PolyLatexIsZero_( coeff )
  try
    tf = isAlways( coeff == 0 );
  catch
    tf = isequal( coeff, sym( 0 ) );
  end
end

function tf = PolyLatexIsOneAbs_( coeff )
  try
    tf = isAlways( coeff == 1 ) || isAlways( coeff == -1 );
  catch
    tf = isequal( coeff, sym( 1 ) ) || isequal( coeff, sym( -1 ) );
  end
end

function [tf, abs_coeff_latex] = PolyLatexSignLatex_( coeff )
  coeff_latex = strtrim( string( latex( coeff ) ) );
  tf = startsWith( coeff_latex, "-" );
  if tf
    abs_coeff = simplify( -coeff );
  else
    abs_coeff = coeff;
  end
  abs_coeff_latex = strtrim( string( latex( abs_coeff ) ) );
end
