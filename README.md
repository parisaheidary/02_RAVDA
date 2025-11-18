P_M_computing_param1 and P_M_computing_S0_param1 are the codes used to generate PODs and the reduced model. The first code, P_M_computing_param1, is used when S0 is known, and the unknowns are parameters. The second code, P_M_computing_S0_param1, is used when S0 and the parameters are coupled.

Reduced_adjoint_script_S0_params includes the minimization functions:

Reduced_Adjoint_fun_S0_params: Calculates J and grad J for S0 using the best guess of parameters (first step) and then uses fminunc to minimize the cost function.
Reduced_Adjoint_fun_S0_params2: Calculates J and grad J for the parameters using the optimal S0 calculated in the first step, then minimizes it using fminunc. This loop repeats until convergence.

Reduced_adjoint_script_param is used when S0 is known and only the parameters are estimated. The function that calculates J and grad J in this case is Reduced_Adjoint_fun_param.

Hydrus_Run_param runs the Hydrus model. It first manipulates the input data and then calls Hydrus to run.

reading_flx_Hy reads the fluxes from ASCII files created after running Hydrus (Hydrus results).

