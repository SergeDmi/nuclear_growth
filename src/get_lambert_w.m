function M = get_lambert_w(x,accuracy)
% Lambert W function is the inverse of x = M*exp(M).

% desired accuracy
if nargin<2
  accuracy=1.e-9;
end

% Initialization
M = ones(size(x));
m = inf*M;

% from Haley
while any(abs(M - m)./abs(M) > accuracy)
   m = M;
   eM = exp(M);
   b = M.*eM - x;
   M = M - b./((eM.*(M+1)-(M+2).*b./(2*M+2)));
end

end
