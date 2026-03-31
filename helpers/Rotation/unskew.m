function v = unskew( vx )
%   unskew takes a 3x3xn matrix of skew symmetric matrix and returns the
%   vector that parameterizes it
%
%   e.g. v = [x y z]', if vx = [ 0 -z  y]
%                              [ z  0 -x]
%                              [-y  x  0]
%
%---------------------------------INPUTS-----------------------------------
%
%   vx:
%       3x3xn array of skew symmetric matrices
%
%--------------------------------OUTPUTS-----------------------------------
%
%   v:
%       3xn array of column vectors
%
%--------------------------------------------------------------------------

%% UNSKEW

%verify proper inputs
[vr,vc,n] = size(vx);
if vr ~= 3 || vc ~= 3
    error('vx must be 3x3xn')
end

%for each vector
v = zeros(3,n);
for k = 1:n
    v(1,k) = vx(3,2,k);
    v(2,k) = vx(1,3,k);
    v(3,k) = vx(2,1,k);       
end



end
