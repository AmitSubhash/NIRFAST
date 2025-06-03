function [pout, vout, viout, facevout, faceiout]  = select3d(obj)

% [pout, vout, viout, facevout, faceiout]  = select3d(obj)
%
% Copyright (c) 2009, The MathWorks, Inc.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the The MathWorks, Inc. nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.


%SELECT3D(H) Determines the selected point in 3-D data space.
%  P = SELECT3D determines the point, P, in data space corresponding 
%  to the current selection position. P is a point on the first 
%  patch or surface face intersected along the selection ray. If no 
%  face is encountered along the selection ray, P returns empty.
%
%  P = SELECT3D(H) constrains selection to graphics handle H and,
%  if applicable, any of its children. H can be a figure, axes, 
%  patch, or surface object.
%
%  [P V] = SELECT3D(...), V is the closest face or line vertex 
%  selected based on the figure's current object. 
%
%  [P V VI] = SELECT3D(...), VI is the index into the object's 
%  x,y,zdata properties corresponding to V, the closest face vertex 
%  selected.
%
%  [P V VI FACEV] = SELECT3D(...), FACE is an array of vertices 
%  corresponding to the face polygon containing P and V. 
%  
%  [P V VI FACEV FACEI] = SELECT3D(...), FACEI is the row index into 
%  the object's face array corresponding to FACE. For patch 
%  objects, the face array can be obtained by doing 
%  get(mypatch,'faces'). For surface objects, the face array 
%  can be obtained from the output of SURF2PATCH (see 
%  SURF2PATCH for more information).
%
%  RESTRICTIONS:
%  SELECT3D supports surface, patch, or line object primitives. For surface 
%  and patches, the algorithm assumes non-self-intersecting planar faces. 
%  For line objects, the algorithm always returns P as empty, and V will
%  be the closest vertex relative to the selection point. 
%
%  Example:
%
%  h = surf(peaks);
%  zoom(10);
%  disp('Click anywhere on the surface, then hit return')
%  pause
%  [p v vi face facei] = select3d;
%  marker1 = line('xdata',p(1),'ydata',p(2),'zdata',p(3),'marker','o',...
%                 'erasemode','xor','markerfacecolor','k');
%  marker2 = line('xdata',v(1),'ydata',v(2),'zdata',v(3),'marker','o',...
%                 'erasemode','xor','markerfacecolor','k');
%  marker2 = line('erasemode','xor','xdata',face(1,:),'ydata',face(2,:),...
%                 'zdata',face(3,:),'linewidth',10);
%  disp(sprintf('\nYou clicked at\nX: %.2f\nY: %.2f\nZ: %.2f',p(1),p(2),p(3)'))
%  disp(sprintf('\nThe nearest vertex is\nX: %.2f\nY: %.2f\nZ: %.2f',v(1),v(2),v(3)'))
% 
%  Version 1.3 11-11-04
%  Copyright Joe Conti 2004 
%  Send comments to jconti@mathworks.com
%
%  See also GINPUT, GCO.

% Output variables
pout = [];
vout = [];
viout = [];
facevout = [];
faceiout = [];

% other variables
ERRMSG = 'Input argument must be a valid graphics handle';
isline = logical(0);
isperspective = logical(0);

% Parse input arguments
if nargin<1
   obj = gco;
end

if isempty(obj) | ~ishandle(obj) | length(obj)~=1
    error(ERRMSG);
end

% if obj is a figure
if strcmp(get(obj,'type'),'figure')
    fig = obj;
    ax = get(fig,'currentobject');
    currobj = get(fig,'currentobject');
    
    % bail out if not a child of the axes
    if ~strcmp(get(get(currobj,'parent'),'type'),'axes')
        return;
    end
    
% if obj is an axes
elseif strcmp(get(obj,'type'),'axes')
    ax = obj;
    fig = get(ax,'parent');
    currobj = get(fig,'currentobject');
    currax = get(currobj,'parent');
    
    % Bail out if current object is under an unspecified axes
    if ~isequal(ax,currax)
        return;
    end

% if obj is child of axes    
elseif strcmp(get(get(obj,'parent'),'type'),'axes')
    currobj = obj;
    ax = get(obj,'parent');
    fig = get(ax,'parent');
    
% Bail out
else 
    return
end

axchild = currobj;
obj_type = get(axchild,'type');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get vertex, face, and current point data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp = get(ax,'currentpoint')';

% If surface object
if strcmp(obj_type,'surface')
    % Get surface face and vertices
    fv = surf2patch(axchild);
    vert = fv.vertices;
    faces = fv.faces;
    
% If patch object
elseif strcmp(obj_type,'patch')
    vert = get(axchild,'vertices');
    faces = get(axchild,'faces');
    
% If line object     
elseif strcmp(obj_type,'line')
    xdata = get(axchild,'xdata');
    ydata = get(axchild,'ydata');
    zdata = get(axchild,'zdata');
    vert = [xdata', ydata',zdata'];
    faces = []; 
    isline = logical(1);

% Ignore all other objects
else     
   return;
end
    
% Add z if empty
if size(vert,2)==2
   vert(:,3) = zeros(size(vert(:,2)));
   if isline
       zdata = vert(:,3);
   end
end

% NaN and Inf check 
nan_inf_test1 = isnan(faces) | isinf(faces);
nan_inf_test2 = isnan(vert) | isinf(vert);
if any(nan_inf_test1(:)) | any(nan_inf_test2(:))
    warning(sprintf('%s does not support NaNs or Infs in face/vertex data.',mfilename));
end

% For debugging
% if 0    
%     ax1 = getappdata(ax,'testselect3d');
%     if isempty(ax1) | ~ishandle(ax1)
%         fig = figure; 
%         ax1 = axes;
%         axis(ax1,'equal');
%         setappdata(ax,'testselect3d',ax1);
%     end
%     cla(ax1);
%     patch('parent',ax1,'faces',faces,'vertices',xvert','facecolor','none','edgecolor','k'); 
%     line('parent',ax1,'xdata',xcp(1,2),'ydata',xcp(2,2),'zdata',0,'marker','o','markerfacecolor','r','erasemode','xor');
% end

% Determine the selection ray
cp1 = cp(:,1);      % point on near plane
cp2 = cp(:,2);      % point on far plane
ray_dir = cp2 - cp1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alternative algorithm that does not rely on hidden properties %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isline
    % Find closest vertex to the selection ray in 3-D
    ray_dir = ray_dir./norm(ray_dir);
    dif = vert - cp1;
    t = dif*ray_dir';
    proj = cp1 + t.*ray_dir;
    d = sqrt(sum((vert-proj).^2,2));
    [~,i] = min(d);
    vout = vert(i,:);
    viout = i;
    return
end

% Initialize intersection search
best_t = inf;
faceiout = [];
viout = [];
vout = [];
pout = [];

ray_dir = ray_dir./norm(ray_dir);

for fi = 1:size(faces,1)
    f = faces(fi,~isnan(faces(fi,:)));
    if numel(f) < 3
        continue
    end
    tris = [f(1) f(2) f(3)];
    [hit, tval] = local_rayTri(cp1, ray_dir, vert(tris,:));
    if ~hit && numel(f)==4
        tris = [f(1) f(3) f(4)];
        [hit, tval] = local_rayTri(cp1, ray_dir, vert(tris,:));
    end
    if hit && tval < best_t && tval > 0
        best_t = tval;
        faceiout = fi;
        pout = cp1 + ray_dir.*tval;
        facevout = vert(f,:);
        [~,idx] = min(sum((facevout - pout).^2,2));
        viout = f(idx);
        vout = vert(viout,:);
    end
end

if isempty(pout)
    return
end
end

function [hit,t] = local_rayTri(orig,dir,tri)
% Ray/triangle intersection using Moller-Trumbore algorithm
eps_val = 1e-9;
v1 = tri(1,:); v2 = tri(2,:); v3 = tri(3,:);
edge1 = v2 - v1;
edge2 = v3 - v1;
h = cross(dir, edge2);
a = dot(edge1,h);
if abs(a) < eps_val
    hit = false; t = Inf; return
end
f = 1.0/a;
s = orig - v1;
u = f*dot(s,h);
if u < 0 || u > 1
    hit = false; t = Inf; return
end
q = cross(s, edge1);
v = f*dot(dir,q);
if v < 0 || (u + v) > 1
    hit = false; t = Inf; return
end
t = f*dot(edge2,q);
if t <= eps_val
    hit = false; t = Inf; return
end
hit = true;
end

