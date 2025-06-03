function pathloc = get_pathloc(loc)

% pathloc = get_pathloc(loc)
%
% determined whether the given location is a variable or path
%
% loc is the location given
% pathloc is the resulting location (in quotes if a path)

if ~contains(loc,'/') && ~contains(loc,'\') && ~contains(loc,'.')
    pathloc = loc;
else
    pathloc = strcat('''',loc,'''');
end