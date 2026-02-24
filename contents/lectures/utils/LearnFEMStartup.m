function LearnFEMStartup()

%% Set all interpreters from tex to latex.
%#ok<*CLALL>
list_factory = fieldnames( get( groot, "factory" ) );
index_interpreter = find( contains( list_factory, "Interpreter" ) );
for ii = 1:length( index_interpreter )
  default_name = strrep( list_factory{index_interpreter(ii)}, "factory", "default" );
  set( groot, default_name, "latex" );
end

%% Set all X/Y/Z LimitMethod to tight
%#ok<*CLALL>
list_factory = fieldnames( get( groot, "factory" ) );
index_interpreter = find( contains( list_factory, "LimitMethod" ) );
for ii = 1:length( index_interpreter )
  default_name = strrep( list_factory{index_interpreter(ii)}, "factory", "default" );
  set( groot, default_name, "tight" );
end

%% Set Box to On
set( groot, "defaultAxesBox", "on" );

%% Clear all startup variables
clear

end

