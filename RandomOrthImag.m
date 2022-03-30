
function [OrthMatrix] = RandomOrthImag(n)
% Schmidt Orthogonalization
% 用施密特正交化输出随机的正交基 大小为 n*n

    in = rand(n)*1i;

    % Init OrthMatrix
    OrthMatrix = zeros(size(in));
    OrthMatrix(:,1) = in(:,1);

    % Orthogonalization
    for k=2:n
        for t=1:k-1
            OrthMatrix(:,k) = OrthMatrix(:,k) - dot(OrthMatrix(:,t), in(:,k))/dot(OrthMatrix(:,t),OrthMatrix(:,t))*OrthMatrix(:,t);
        end
        OrthMatrix(:,k) = OrthMatrix(:,k) + in(:,k);
    end

    % Normalization
    for k=1:n
        OrthMatrix(:,k) = OrthMatrix(:,k)/norm(OrthMatrix(:,k), 'fro');
    end
end