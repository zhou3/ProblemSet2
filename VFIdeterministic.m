close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
T_mat = [0.977, 1-0.977; 1 - 0.926, 0.926];
A = [1.1, 0.678];

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 11; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
cons(:,:,1) = A(1) * k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 
cons(:,:,2) = A(2) * k_mat .^ alpha + (1 - delta) * k_mat - k_mat';

ret = cons .^ (1 - sigma) / (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret(cons < 0) = -Inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(2, num_k);
while dis > tol
    % an alternative, more direct way to compute the value array:
    value_mat_alt = ret + beta * ...
        repmat(permute((T_mat * v_guess), [3 2 1]), [num_k 1 1]);
    
    % compute the utility value for all possible combinations of k and k':
    value_mat(:,:,1) = ret(:,:,1) + beta * ( ...
        T_mat(1,1) * repmat(v_guess(1,:), [num_k 1]) + ...
        T_mat(1,2) * % finish from here!
    
    
    % find the optimal k' for every k:
    [vfn, pol_indx] = max(value_mat, [], 2);
    vfn = vfn';
    
    % what is the distance between current guess and value function
    dis = max(abs(vfn - v_guess));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
end

g = k(pol_indx); % policy function

plot(k,vfn)
figure
plot(k,g)


