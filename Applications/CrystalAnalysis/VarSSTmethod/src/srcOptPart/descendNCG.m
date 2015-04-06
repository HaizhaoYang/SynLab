function [x,e,de] = descendNCG( f, x, maxNumSteps, display )
  stepsize = 1;
  [e,de] = f( x );
  p = 0 * de;
  normSqrOld = 1;
  %% NCG iteration
  for ii = 1:maxNumSteps
    if display
      fprintf('step %d, stepsize %f, energy %f, gradient-norm %e\n',ii-1,stepsize,e,norm(de(:),2))
    end
    %% compute search direction
    normSqr = sum( de(:).^2 );
    p = -de + normSqr / normSqrOld * p;
    expDec = sum( de(:) .* p(:) );
    if expDec >= 0
      p = -de;
      expDec = sum( de(:) .* p(:) );
    end
    %% do descent step
    stepsize = min( stepsize, - e / expDec );
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
    de = trialDE;
    normSqrOld = normSqr;
  end
end
