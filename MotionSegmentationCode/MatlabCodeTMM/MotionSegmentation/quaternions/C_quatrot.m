function Y = C_quatrot(X,Q)
% Y = quatrot(X,Q)
% Rotate vectors in 3xN matrix X by quaternions in Q.
% If Q is 4xN, "point-wise" rotation will be performed.
% If Q is 4x1, all vectors in X will be rotated by Q.
%
% Remark:  To concatenate rotations in quaternion representation, multiply them from left to right, 
%          just like rotation matrices.
% Example: First rotating by q1, then by q2 is equivalent to rotating by quatmult(q2,q1).

if (size(X,1)~=3)
    error('X has to be 3xN.');
end
if (size(Q,1)~=4)
    error('Q has to have 4 rows (quaternions!).');
end

N = size(X,2);

if (size(Q,2)>1) && (size(X,2)>1 && size(X,2)~=size(Q,2))
    error('Size mismatch between X and Q!');
end
padder = [zeros(1,N);X];
if(size(padder,2)==1)
    padder = repmat(padder,1,size(Q,2));
end

Y = quatmult(quatmult(Q,padder),quatinv(Q));
Y = Y(2:4,:);
