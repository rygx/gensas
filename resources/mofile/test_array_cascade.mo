// A cascade of first-order lags written with an array and a for-equation.
// x[1] = exp(-t), x[2] = t exp(-t), x[3] = t^2/2 exp(-t), ...
model Cascade
  parameter Integer n = 5;
  Real x[n](start = {1, 0, 0, 0, 0});
equation
  der(x[1]) = -x[1];
  for i in 2:n loop
    der(x[i]) = x[i-1] - x[i];
  end for;
end Cascade;
