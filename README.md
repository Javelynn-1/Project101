# Project101 MATLAB R2024b False Position
# Thanawich Adisaisakda 6610110121
# Copy below this
-------------------------------------------------
%% Exercise 5.10 - False Position Method

% Define the function whose root we want to approximate
Q = 20;
g = 9.81;
B = @(y) 3 + y; % Function B
Ac = @(y) 3*y + y.^2/2; % Function Ac
func = @(y) (1 - Q^2 * B(y) ./ (g * Ac(y).^3)); % Function to find the root

% Define the upper and lower guesses
xu = 2.5;
xl = 0.5;

% Define stopping criterion and max iterations
es = 1; % Desired relative error (1% tolerance)
n_iter = 23; % Set maximum number of iterations to 23

%% Numerically - False Position method
[root_2, f_2, e_2, n_2, results] = false_position(func, xl, xu, es, n_iter);

% Display results for False Position method
fprintf('False Position method: \n');
fprintf('\t root: %.4f\n', root_2);
fprintf('\t function value at root: %.4f\n', f_2);
fprintf('\t approximate relative error: %.4f\n', e_2);
fprintf('\t number of iterations: %d\n\n', n_2);

% Display results for each iteration
disp('Iteration Results (root, f(root), relative error):');
disp('-------------------------------------------');
disp('Iteration\tRoot\t\tf(Root)\t\tError');
for i = 1:n_2
    fprintf('%d\t\t%.4f\t\t%.4f\t\t%.4f\n', results(i,1), results(i,2), results(i,3), results(i,4));
end

%% Graphically
x = linspace(0.5, 2.5, 100); % x values for plotting
figure;
plot(x, func(x), 'r', 'LineWidth', 2); % Plot func(x)
hold on;
plot(x, zeros(size(x)), 'b--', 'LineWidth', 1.5); % Plot y = 0
scatter(results(:,2), zeros(size(results(:,2))), 'ko', 'filled'); % Plot iteration points
xlabel('y');
ylabel('f(y)');
title('False Position Method - Root Finding');
grid on;
hold off;

%% False Position method implementation
function [root, f_root, ea, iter, all_results] = false_position(func, xl, xu, es, max_iter)
    iter = 0;
    ea = 100; % initial error
    all_results = []; % To store results for each iteration
    xr_old = xl; % Initialize previous root

    while iter < max_iter && ea > es
        % Calculate the root using False Position formula
        xr = xu - (func(xu) * (xl - xu)) / (func(xl) - func(xu));
        f_xr = func(xr);
        
        % Calculate approximate relative error
        if iter > 0
            ea = abs((xr - xr_old) / xr) * 100;
        end
        
        % Store previous root
        xr_old = xr;
        
        % Save results for this iteration
        all_results = [all_results; iter, xr, f_xr, ea];
        
        % Update the interval
        if func(xl) * f_xr < 0
            xu = xr; % Root is in the left subinterval
        elseif func(xr) * func(xu) < 0
            xl = xr; % Root is in the right subinterval
        else
            break; % If f(xr) == 0, stop iteration
        end
        
        iter = iter + 1; % increment iteration count
    end
    
    root = xr; % return the root
    f_root = f_xr; % return the function value at the root
end
