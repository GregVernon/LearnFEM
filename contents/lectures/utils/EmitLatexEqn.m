function EmitLatexEqn( expression_string )
% EmitLatexEqn  Emit equation markup from MATLAB for Quarto (%| output: asis).
%
% By default this emits display math using $$...$$ so it renders in both
% HTML (MathJax) and PDF. Use mode="latex" for LaTeX-only raw output.
%
% Usage:
%   EmitLatexEqn( "p(x) = " + string(latex(p)) );
%   EmitLatexEqn( "\\lambda(C) \\approx " + string(latex(v)) );

arguments
  expression_string (1,1) string
end

  disp( "```{=latex}" )
  disp( "\begin{equation}" )
  disp( "  " + expression_string )
  disp( "\end{equation}" )
  disp( "```" )
  disp( "" )

end
