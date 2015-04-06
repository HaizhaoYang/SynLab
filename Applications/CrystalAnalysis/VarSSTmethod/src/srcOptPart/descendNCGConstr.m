function [x,e,de] = descendNCGConstr( f, x0, P, maxNumSteps, display )
% min f(x) s.t. (I-P)(x0-x)=0, where P is an orthogonal projection onto a subspace
% assumes feasible starting point

  stepsize = 1;
  x = x0;
  [e,de] = f( x );
  de = P(de);
  p = 0 * de;
  normSqrOld = 1;
  %% NCG iteration
  for ii = 1:maxNumSteps
    if display == 2
      fprintf('step %d, stepsize %f, energy %f, gradient-norm %e, constraint %e\n',ii-1,stepsize,e,norm(de(:),2),norm(reshape(x-x0-P(x-x0),[],1),2))
    elseif display == 1
      fprintf('step %d, stepsize %f, energy %f, gradient-norm %e\n',ii-1,stepsize,e,norm(de(:),2))
    end
    %% compute search direction
    normSqr = sum( de(:).^2 );
    p = -de + normSqr / normSqrOld * p;
    expDec = sum( de(:) .* p(:) );
    if expDec >= 0
      p = -de;
      expDec = -expDec;
    end
    %% do descent step
    stepsize = 1;
    trialX = x + stepsize * p;
    [trialE,trialDE] = f( trialX );
    %% Armijo linesearch
    while ( trialE > e + .5 * stepsize * expDec ) || isnan(trialE)
      stepsize = stepsize / 2;
      trialX = x + stepsize * p;
      [trialE,trialDE] = f( trialX );
    end
    x = trialX;
    e = trialE;
    de = P(trialDE);
    if mod(ii,10) == 0
      x = x0 + P(x-x0);
    end
    normSqrOld = normSqr;
  end
end
