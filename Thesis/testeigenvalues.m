% Define the system matrix A and control matrix B
A = all_models{1}; % Define your system matrix
B = B_matrix{1}; % Define your control matrix

% Compute the eigenvalues of the closed-loop system
eigenvalues = eig(A - B*controllerGain{1}); % K is the control law
% Define the system matrix A and control matrix B


% Plot the eigenvalues
scatter(real(eigenvalues), imag(eigenvalues), 'o');
hold on; % Add the unit circle without clearing the previous plot
viscircles([0,0], 1, 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'black');
xlabel('Real Part');
ylabel('Imaginary Part');
title('Eigenvalues of the Closed-Loop System');
grid on;
hold off; % Release the hold on the plot
