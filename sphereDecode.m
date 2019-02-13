function [x_hat, e] = sphereDecode(y,H,R,const,q,r)
%%% Sphere decoder
% Inputs:
%   y - nx1 received vector
%   H - nxn channel gain matrix (or estimate)
%   R - decoding radius
%   const - constellation to search over
%   q - optional Q of QR decomposition of H
%   r - optional R of QR decomposition of H
%
% Outputs:
%   x_hat - estimate of x, will be zero if decoding fails
%   e - L2 norm from x_hat to nearest point

n = length(H);

if nargin < 5
    [q,r] = qr(H);
end

state.R = R;
state.r = r;
state.qy = q'*y;
state.vec = [0; 0; 0; 0];
state.depth = n;
state.best_l2 = 1e15;
state.best_vec = [0; 0; 0; 0];
state.const = const;
state = sphereRecurse( state );
state.best_vec;
state.best_l2;

x_hat = state.best_vec;
e = state.best_l2;
if e == 1e15
    x_hat = zeros(1,n);
end

function state = sphereRecurse( state )
    ii = state.depth;
    for jj = 1:length(state.const)
        state.vec(ii) = state.const(jj);
        l2_part = (state.qy(ii) - state.r(ii,:)*state.vec).^2;
        
        %Only recurse if we are withing the decoding bound
        if l2_part < state.R 
            if ii > 1
                %If we are not yet at the root, then recurse
                state.depth = ii - 1;
                state = sphereRecurse( state );
            else
                %If we are at a root then see if we've found the best point
                if l2_part < state.best_l2
                    state.best_vec = state.vec;
                    state.best_l2 = l2_part;
                end
            end
        end
    end
end

end