
function h = addcircle(h, c, radius, weight)

    [hr, hc] = size(h);
    
    if nargin == 3
	weight = 1;
    end
    
    % c and radius must be integers
    if any(c-fix(c))
	error('Circle centre must be in integer coordinates');
    end
    
    if radius-fix(radius)
	error('Radius must be an integer');
    end
    
    x = 0:fix(radius/sqrt(2));
    costheta = sqrt(1 - (x.^2 / radius^2));
    y = round(radius*costheta);
    
    % Now fill in the 8-way symmetric points on a circle given coords 
    % [px py] of a point on the circle.
    
    px = c(2) + [x  y  y  x -x -y -y -x];
    py = c(1) + [y  x -x -y -y -x  x  y];

    % Cull points that are outside limits
    validx = px>=1 & px<=hr;
    validy = py>=1 & py<=hc;    
    valid = find(validx & validy);

    px = px(valid);
    py = py(valid);
    
    ind = px+(py-1)*hr;
    h(ind) = h(ind) + weight;
