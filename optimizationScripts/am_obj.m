function [r, g] = am_obj(x, varargin)

if ~isempty(varargin)
    r = varargin{1};
else
    r = x(end);
end
g = zeros(size(x));
g(end) = 1;

end