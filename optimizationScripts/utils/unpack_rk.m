function [A, b, r] = unpack_rk(x, s, mthType)
% DIRK-G : x = [A(:), b(:) r]

xx = x;

if strcmpi(mthType, 'erk')
    % explicit method
    an = 0.5*s*(s-1);
    
    xa = xx(1:an);
    A = zeros(s);
    A(tril(true(s),-1)) = xa;
    
elseif strcmpi(mthType, 'dirk')
    % implicit method
    an = 0.5*s*(s+1);
    
    xa = xx(1:an);
    A = zeros(s);
    A(tril(true(s),0)) = xa;
end

try
    % unpack b
    b = xx(an+1:an+s);
    r = nan;
    
catch err
    keyboard
end

end
