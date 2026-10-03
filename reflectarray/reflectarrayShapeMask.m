function keep = reflectarrayShapeMask(X,Y,shape)
%REFLECTARRAYSHAPEMASK Clip lattice sites before skew or manual movement.
%   X and Y are the centred, orthogonal design-lattice coordinates.

halfW = max(abs(X(:)));
halfH = max(abs(Y(:)));
if halfW == 0 || halfH == 0
    keep = true(size(X));
    return;
end
switch char(shape)
    case {'Custom','Rectangle'}
        keep = true(size(X));
    case 'Diamond'
        keep = abs(X)/halfW + abs(Y)/halfH <= 1+1e-9;
    case 'Circle'
        r = min(halfW,halfH);
        keep = X.^2+Y.^2 <= r^2*(1+1e-9);
    case 'Ellipse'
        % Deliberately wider than Circle even on a square lattice.
        rx = max(1.15*halfW,eps);
        ry = max(min(halfH,0.75*halfW),eps);
        keep = (X/rx).^2+(Y/ry).^2 <= 1+1e-9;
    case 'Hexagon'
        keep = polygonMask(X,Y,6,min(halfW,halfH));
    case 'Octagon'
        keep = polygonMask(X,Y,8,min(halfW,halfH));
    otherwise
        error('Reflectarray:Shape','Unknown surface shape.');
end
% Very small even grids may have no centre lattice site inside a mask.
if ~any(keep(:)), keep = true(size(X)); end
end

function keep = polygonMask(X,Y,n,apothem)
keep = true(size(X));
for side = 0:n-1
    angle = side*360/n;
    keep = keep & X*cosd(angle)+Y*sind(angle) <= apothem*(1+1e-9);
end
end
