maxNumCompThreads(1);

%% == Select Problem ================================================================================================ %%
problem_index = 1; % 1 - ADR2D, 2 - Burgers1d

switch problem_index
    case 1
        problem = ADR2d();
        num_timesteps = 250;
    case 2
        problem = Burgers1d();
        num_timesteps = 500;
end

%% == Solvers ======================================================================================================= %%
record_stats = true;
linear_solver_type  = 1; % 1 - mldivide, 2 - gmres

% -- Linear Solver -----------------------------------------------------------------------------------------------------
switch linear_solver_type
    case 1
        ls_constructor  = @MLDV;
        ls_options = struct(                ...
            'record_stats', record_stats    ...
            );
    case 2
    	ls_constructor  = @GMRESMAT;
        ls_options = struct(                    ...
            'max_iterations',   100,            ...
            'tolerance',        1e-12,          ...
            'record_stats',     record_stats    ...
            );
end

% -- Nonlinear Solver --------------------------------------------------------------------------------------------------
record_stats = true;
ls_handle  = @() ls_constructor(ls_options);
nls_handle = @() Newton(struct('linear_solver', ls_handle(), 'record_stats', record_stats));


%% == Start parallel Pool =========================================================================================== %%
spmd_burnin(problem);

%% == BDF Integrator ================================================================================================ %%
bdf_order = 4;
starting_integrator = epirk4s3Const();
bdf = BDFConst( ...
    struct( ...
        'order',                bdf_order, ...
        'nonlinear_solver',     nls_handle(), ...
        'starting_integrator',  starting_integrator, ... 
        'record_stats',         true ...
    ) ...
);
bdf.num_timesteps = num_timesteps;

%% == BBDF Integrators ============================================================================================== %%
bbdf_order = 4;

% -- serial bbdf -------------------------------------------------------------------------------------------------------
bbdf = BBDFConst( ...
    	struct( ...
            'z',                    1i * linspace(-1, 1, bbdf_order), ...
            'alpha',                1 / 2, ...
            'nonlinear_solver',     nls_handle(), ...
            'starting_integrator',  starting_integrator, ...
            'record_stats',         true ...
        ) ...
    );
bbdf.num_timesteps = num_timesteps;

% -- parallel bbdf (requires matlab parallel toolkit) ------------------------------------------------------------------
bbdf_spmd = BBDFConst_spmd( ...
    	struct( ...
            'z',                    1i * linspace(-1, 1, bbdf_order), ...
            'alpha',                1 / 2, ...
            'nonlinear_solver',     nls_handle(), ...
            'starting_integrator',  starting_integrator, ...
            'record_stats',         true ...
        ) ...
    );
bbdf_spmd.num_timesteps = num_timesteps;

%% == AM Integrator ================================================================================================= %%
am_order = 2;
starting_integrator = epirk4s3Const();
am = AMConst( ...
        struct( ...
            'order', am_order, ...
            'nonlinear_solver', nls_handle(), ...
            'starting_integrator', starting_integrator, ...
            'record_stats', true ...
        ) ...
    );
am.num_timesteps = num_timesteps;

%% == BAM Integrators =============================================================================================== %%
bam_order = 4;

% -- serial bam --------------------------------------------------------------------------------------------------------
bam = BAMConst( ...
    	struct( ...
            'z',                    1i * linspace(-1, 1, bam_order), ...
            'b',                    1i * linspace(-1, 1, bam_order), ...
            'alpha',                1 / 2, ...
            'nonlinear_solver',     nls_handle(), ...
            'starting_integrator',  starting_integrator, ...
            'record_stats',         true ...
        ) ...
    );
bam.num_timesteps = num_timesteps;

% -- parallel bam (requires matlab parallel toolkit) -------------------------------------------------------------------
bam_spmd = BAMConst_spmd( ...
    	struct( ...
            'z',                    1i * linspace(-1, 1, bam_order), ...
            'b',                    1i * linspace(-1, 1, bam_order), ...
            'alpha',                1 / 2, ...
            'nonlinear_solver',     nls_handle(), ...
            'starting_integrator',  starting_integrator, ...
            'record_stats',         true ...
        ) ...
    );
bam_spmd.num_timesteps = num_timesteps;


%% == Example Solve ================================================================================================= %%

y_bdf = bdf.solve(problem);
fprintf('error bdf: %f\ttime bdf: %f\n', max(abs(y_bdf - problem.reference_solution)), bdf.step_stats.total_seconds)

y_am = am.solve(problem);
fprintf('error am: %f\ttime am: %f\n', max(abs(y_am - problem.reference_solution)), am.step_stats.total_seconds)

y_bbdf = bbdf.solve(problem);
fprintf('error bbdf: %f\ttime bbdf: %f\n', max(abs(y_bbdf - problem.reference_solution)), bbdf.step_stats.total_seconds)

y_bam = bam.solve(problem);
fprintf('error bam: %f\ttime bam: %f\n', max(abs(y_bam - problem.reference_solution)), bam.step_stats.total_seconds)

y_bbdf_spmd = bbdf_spmd.solve(problem);
fprintf('error bbdf_spmd: %f\ttime bbdf_spmd: %f\n', max(abs(y_bbdf_spmd - problem.reference_solution)), bbdf_spmd.step_stats.total_seconds)

y_bam_spmd = bam_spmd.solve(problem);
fprintf('error bam_spmd: %f\ttime bam_spmd: %f\n', max(abs(y_bam_spmd - problem.reference_solution)), bam_spmd.step_stats.total_seconds)