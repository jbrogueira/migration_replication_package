function [Vsol, info] = anderson_solve(Vstart, opts, np)
%% Anderson-accelerated fixed-point solve of the transition path: find V with
%% V = Phi(V) (sweep_phi). Each iteration costs ONE sweep. Window m residual
%% differences; damped mixing beta_mix on the accelerated update. Suited to
%% the open-budget system (mild spectrum); not for the closed-budget map.
  if nargin < 3; np = struct(); end
  d = @(f, v) merge(np, f, v);
  tol_F = d('tol_F', 5e-5); maxit = d('maxit', 400);
  mwin = d('m', 15); bmix = d('beta_mix', 0.5);
  tol_win = d('tol_win', 0); win_slots = d('win_slots', []);

  tracks = {'V_eu15_h','V_nms8_h','V_oth_h','V_eu15_l','V_nms8_l','V_oth_l', ...
            'V_eu15_h_std','V_nms8_h_std','V_oth_h_std', ...
            'V_eu15_h_pref','V_nms8_h_pref','V_oth_h_pref'};
  time = size(Vstart.(tracks{1}), 2); nr = size(Vstart.(tracks{1}), 1); ns = time-2;
  n = numel(tracks)*nr*ns;
  winmask = [];
  if ~isempty(win_slots)
      base = false(nr, ns); base(:, (win_slots(1)-1):(win_slots(2)-1)) = true;
      winmask = repmat(base(:), numel(tracks), 1);
  end

  stack = @(S) cell2mat(cellfun(@(t) reshape(S.(t)(:,2:time-1), [], 1), tracks', 'UniformOutput', false));
  function S = unstack(v)
      S = Vstart; p = 0;
      for q = 1:numel(tracks)
          S.(tracks{q})(:,2:time-1) = reshape(v(p+1:p+nr*ns), nr, ns); p = p + nr*ns;
      end
  end

  x = stack(Vstart);
  X = []; G = [];   % histories of x and g(x) = Phi(x)
  hist = [];
  for it = 1:maxit
      [P, dg] = sweep_phi(unstack(x), opts);
      g = stack(P);
      f = g - x;                       % residual
      nF = max(abs(f));
      hist(end+1) = nF;  %#ok<AGROW>
      nFw = nF; if ~isempty(winmask); nFw = max(abs(f(winmask))); end
      if mod(it, 10) == 1 || nF < tol_F || (tol_win>0 && nFw < tol_win)
          printf('anderson %3d: |F|_inf = %.6e  window |F| = %.6e\n', it-1, nF, nFw); fflush(stdout);
      end
      if nF < tol_F; break; end
      if tol_win > 0 && nFw < tol_win; printf('window tolerance met.\n'); break; end
      X = [X, x]; G = [G, g];          %#ok<AGROW>
      if size(X, 2) > mwin; X(:,1) = []; G(:,1) = []; end
      k = size(X, 2);
      if k == 1
          x = x + bmix*f;
      else
          dF = diff(G - X, 1, 2);      % residual differences
          gamma = dF \ f;              % least squares
          xbar = X(:,end) - (X(:,2:end) - X(:,1:end-1))*gamma;
          gbar = G(:,end) - (G(:,2:end) - G(:,1:end-1))*gamma;
          x = (1-bmix)*xbar + bmix*gbar;
      end
  end
  Vsol = unstack(x);
  info = struct('F_final', hist(end), 'F_hist', hist, 'iters', numel(hist), 'dg', dg);
end
function v = merge(s, f, v0)
  if isfield(s, f); v = s.(f); else; v = v0; end
end
