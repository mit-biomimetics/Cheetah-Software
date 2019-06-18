% Import functionality for Octave
% Place in Octave's path, but never in MATLAB's path as it conflicts with,
% and is inferior to, MATLAB's built-in IMPORT function.
function import(varargin)
 % Import a package. Only entire packages can be imported currently.
 error(nargchk(1, inf, nargin, 'struct'));
 % Import the packages one-by-one
 for i=1:nargin
   import1(varargin{i});
 end
end

function import1(pkgname)
  % Split up name into parts
  pkgname_parts = strsplit(pkgname, '.');
  % Fallback to system import command if not entire package (i.e. failure)
  if length(pkgname_parts)~= 2 || ~strcmp(pkgname_parts{end}, '*')
    error('Only the syntax ''import package_name.*'' is currently supported')
  end
  % Get path for package
 pkgpath = locatepkg(pkgname_parts{1});
  % Add to path
  addpath(pkgpath);
end

function pkgpath = locatepkg(pkgname)
 pathdirs = strsplit(path, pathsep);
 for iPath=1:length(pathdirs)
   pkgpath = [pathdirs{iPath} filesep '+' pkgname];
   if exist(pkgpath, 'dir')
     return;
   end
 end
 error('Package ''%s'' cannot be located in the path', pkgname);
end
