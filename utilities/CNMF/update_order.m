function [O, lo] = update_order(A)

K = size(A, 2); % Number of components
F = (A' * A > 0);  % Create an overlap matrix: components with non-zero overlap
F(1:K+1:K^2) = 0;  % Remove diagonal elements (no self-overlap)
rem_ind = 1:K;      % Initialize remaining indices

% Initialize O and lo
O = {};
lo = [];

if nnz(F) == 0  % Check if there are no overlaps (F is all zeros)
    % No overlaps: return a default order (1:K)
    O{1} = 1:K;
    lo = K;
else
    dp = 0; % Iteration counter
    % While there are remaining indices, update the order
    while ~isempty(rem_ind)
        dp = dp + 1;
        % Select indices with the smallest overlap, using vertex cover
        L = sort(app_vertex_cover(F(rem_ind, rem_ind)), 'ascend');
        ord_ind = setdiff(rem_ind, rem_ind(L)); % Find remaining indices after removing those with overlap
        O{dp} = ord_ind; % Store the current order of components
        lo(dp) = length(ord_ind); % Store the number of non-overlapping components
        rem_ind = rem_ind(L); % Update remaining indices
    end

    % Reverse the order if it's not empty
    if ~isempty(O)
        O = fliplr(O);
        lo = fliplr(lo);
    end
end
end

function L = app_vertex_cover(A)
L = [];
while ~isempty(find(A, 1))
    i = randsample(find(A),1);
    [u,~] = ind2sub(size(A),i);
    L = [L,u];
    A(u,:) = 0; % edges from u
    A(:,u) = 0; % edges to u
end
end


%
% function [O,lo] = update_order(A)
%
%     K = size(A,2);
%     F = (A'*A>0);       % find overlapping components
%     F(1:K+1:K^2) = 0;   % remove diagonal elements
%     rem_ind = 1:K;      % remaining indices
%
%     dp = 0;
%     while ~isempty(rem_ind);
%         dp = dp+1;
%         L = sort(app_vertex_cover(F(rem_ind,rem_ind)),'ascend');
%         ord_ind = setdiff(rem_ind,rem_ind(L));
%         O{dp} = ord_ind;
%         lo(dp) = length(ord_ind);
%         rem_ind = rem_ind(L);
%     end
%     O = fliplr(O);
%     lo = fliplr(lo);
%
%     function L = app_vertex_cover(A)
%         L = [];
%         while ~isempty(find(A, 1))
%             i = randsample(find(A),1);
%             [u,~] = ind2sub(size(A),i);
%             L = [L,u];
%             A(u,:) = 0; % edges from u
%             A(:,u) = 0; % edges to u
%         end
%     end
% end