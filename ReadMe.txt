*** This file explain the code structure in this folder ***

*** The code is developed to implement ALM to solve a general convex problem ***

*** Algorithm Overview ***

The algorithm is implemented under double floating point data and arithmetics (default setting in the Matlab 2015a on this machine)

Given problem: min f(x) s.t. Ax = b, x \in X

Build Augmented Lagrangian function: 

L_\rho(x, \lambda) = f(x) + <\lambda, Ax-b> + \rho/2 * \|Ax - b\|^2

Then given initial multiplier \lambda_0, iteratively implement following steps:

	Step 1: x_k = argmin_{x \in X} L_\rho(x, \lambda_k)

	Step 2: \lambda_{k+1} = \lambda_{k} + \rho*(Ax_k - b)

Until \|Ax_k - b\| is smaller than some user-specific threshold \epsilon

*** Code structure ***

* init_problem.m *
— specify the problem type
- create the problem data in a data struct in the workspace
- solve the problem using standard methods, according to problem type

* init_param_ALM.m *
- specify the algorithm parameter
- assign variable space
- specify how to solve the inner problem
- create a param and a var structs in the workspace

* inner_closeform.m *
- solve inner problem using closed form solution
- check problem specification, if not solvable in closed form, send an alert
- if solvable, output the solution according to problem types

* inner_GPM.m *
- solve the inner problem using GPM (with fixed step size)
- specify how to terminate the problem: 1. by progress at each iteration 2. by specific number of iteration (even goes with 1, still need to assign a max iter number)
- calculate step size according to problem type and data
- assign variable space
- run GPM algorithms
- need to call functions: calc_grad, calc_proj;
- calculate the ergodic average

* outer_update(no .m file) *
- check stopping criteria, for example, \|Ax_k - b\| < \epsilon
- if not stopped, update the multiplier (without projection)

* makeplot.m*






