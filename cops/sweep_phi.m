function [Vout, dg] = sweep_phi(Vin, opts)
%% One UNDAMPED sweep of the Approach 2 transition map (PATH_SOLVER_SPEC.md sec 1):
%% value-track state in -> forward migration/stocks -> temporary equilibrium ->
%% budget clearing -> wedges -> backward value update -> new value-track state out.
%% Vin/Vout: struct with the 12 value tracks (V_eu15_h, V_nms8_h, V_oth_h,
%% V_eu15_l, V_nms8_l, V_oth_l, V_*_h_std, V_*_h_pref), each (J*N or 34) x time.
%% opts fields (all optional): budget_off, pref_sched (logical); tau_p_pref,
%% time_horizon, regain_priced (numeric); diag_dir (string, default 'iter_diag_sweep').
%% dg: diagnostics — Ydiff/Ydiff_full per slot, Ymax, lambda_hat_path, f_std_path,
%% f_l_path, f_pref_path, budget-residual path, plus J, N, time.
%% The model equations live in sweep_engine.m (an exact copy of
%% COPS_shock_europe_skills.m plus two marked insertions). This wrapper only
%% sets the environment, injects the state, runs ONE sweep, and collects Ynew
%% (the pre-damping output) — it implements Phi, not the damped update.

  if nargin < 2; opts = struct(); end
  flag = @(name, on) merge_env(name, on);
  flag('OUTER_MAX', '1');
  flag('BUDGET_OFF',  cond_str(isfield(opts,'budget_off') && opts.budget_off));
  flag('PREF_SCHED',  cond_str(isfield(opts,'pref_sched') && opts.pref_sched));
  %% Baseline-schedule treatment: engine default is INDEXED (Treatment 2, author
  %% decision 2026-07-02). Pass opts.tau_indexed = false to run the unindexed
  %% comparison (schedule fixed in wage units).
  if isfield(opts,'tau_indexed') && ~opts.tau_indexed
      setenv('TAU_INDEXED','0');
  else
      unsetenv('TAU_INDEXED');
  end
  if isfield(opts, 'tau_p_pref');   setenv('TAU_P_PREF',   num2str(opts.tau_p_pref, '%.10g')); else; unsetenv('TAU_P_PREF'); end
  if isfield(opts, 'time_horizon'); setenv('TIME_HORIZON', num2str(opts.time_horizon)); else; unsetenv('TIME_HORIZON'); end
  if isfield(opts, 'diag_dir'); setenv('DIAG_DIR', opts.diag_dir); else; setenv('DIAG_DIR', 'iter_diag_sweep'); end
  if isfield(opts, 'tempeq_tol'); setenv('TEMPEQ_TOL', num2str(opts.tempeq_tol, '%.3e')); else; unsetenv('TEMPEQ_TOL'); end
  if isfield(opts, 'pref_mix'); setenv('PREF_MIX', num2str(opts.pref_mix, '%.10g')); else; unsetenv('PREF_MIX'); end
  if isfield(opts, 'agglom_a'); setenv('AGGLOM_A', num2str(opts.agglom_a, '%.10g')); else; unsetenv('AGGLOM_A'); end
  if isfield(opts, 'tempeq_maxit'); setenv('TEMPEQ_MAXIT', num2str(opts.tempeq_maxit)); else; setenv('TEMPEQ_MAXIT', '20000'); end  %% hang guard default for solver use
  if isfield(opts, 'k_temp_pref'); setenv('K_TEMP_PREF', num2str(opts.k_temp_pref)); else; unsetenv('K_TEMP_PREF'); end
  if isfield(opts, 'target_ctry'); setenv('TARGET_CTRY', num2str(opts.target_ctry)); else; unsetenv('TARGET_CTRY'); end  %% target-country index (default PT=16)
  if isfield(opts, 'regain_priced'); setenv('REGAIN_PRICED', num2str(opts.regain_priced)); else; unsetenv('REGAIN_PRICED'); end  %% price regain of eligibility at the stock rate lambda_r2 (default off)

  INJECT_V   = Vin;  %#ok<NASGU>  % read by sweep_engine insertion 1
  if isfield(opts, 'lambda_path'); INJECT_LAMBDA = opts.lambda_path; end  %#ok<NASGU>
  SWEEP_MODE = 1;    %#ok<NASGU>  % engine raises sweep:done after the loop (insertion 2)
  try
      run('sweep_engine.m');
      error('sweep_phi:noStop', 'sweep_engine ran to completion without the sweep:done stop');
  catch err
      if ~strcmp(err.identifier, 'sweep:done'); rethrow(err); end
  end

  %% Ynew: stacked unsplit (6*J*N x time), already carries the terminal slice = 1.
  %% Split sweep output: Ynew_h_std / Ynew_h_pref (3*J*N x time).
  Vout = struct();
  Vout.V_eu15_h = Ynew(1:J*N, :);
  Vout.V_nms8_h = Ynew(J*N+1:2*J*N, :);
  Vout.V_oth_h  = Ynew(2*J*N+1:3*J*N, :);
  Vout.V_eu15_l = Ynew(3*J*N+1:4*J*N, :);
  Vout.V_nms8_l = Ynew(4*J*N+1:5*J*N, :);
  Vout.V_oth_l  = Ynew(5*J*N+1:6*J*N, :);
  Vout.V_eu15_h_std  = Ynew_h_std(1:J*N, :);
  Vout.V_nms8_h_std  = Ynew_h_std(J*N+1:2*J*N, :);
  Vout.V_oth_h_std   = Ynew_h_std(2*J*N+1:3*J*N, :);
  Vout.V_eu15_h_pref = Ynew_h_pref(1:J*N, :);
  Vout.V_nms8_h_pref = Ynew_h_pref(J*N+1:2*J*N, :);
  Vout.V_oth_h_pref  = Ynew_h_pref(2*J*N+1:3*J*N, :);
  %% Raw passthrough: slot 1 and the terminal slice are returned exactly as the
  %% engine produced them (slot 1 is never written by the backward update and
  %% stays zero; the terminal slice is set to 1 inside the engine). Callers
  %% (regression driver, Newton residual) own the boundary conventions.

  dg = struct('Ymax', Ymax, 'Ydiff', Ydiff, 'Ydiff_full', Ydiff_full, ...
              'lambda_hat_path', lambda_hat_path, 'f_std_path', f_std_path, ...
              'f_l_path', f_l_path, 'f_pref_path', f_pref_path, ...
              'budget_resid', diag_budget_resid, 'J', J, 'N', N, 'time', time, ...
              'denom', diag_denom, 'A_hat', diag_A_hat, 'B_hat', diag_B_hat, ...
              'C_hat', diag_C_hat, 'G_target', diag_G_base_t, ...
              'mult_std_path', mult_std_path, 'mult_pref_path', mult_pref_path, ...
              'mult_l_path', mult_l_path, 'What_h_path', What_h_path, ...
              'What_l_path', What_l_path, ...
              'tau_h_t_path', tau_h_t_path, 'tau_l_t_path', tau_l_t_path);
  %% f_*_path = LEVEL wedges F(t); mult_*_path = the sequential changes F(t)/F(t-1)
  %% actually applied to realwages; What_*_path = counterfactual/baseline wage level
  %% ratios at PT (convention fix 2026-07-02).
end

function s = cond_str(b)
  if b; s = '1'; else; s = ''; end
end

function merge_env(name, val)
  if isempty(val); unsetenv(name); else; setenv(name, val); end
end
