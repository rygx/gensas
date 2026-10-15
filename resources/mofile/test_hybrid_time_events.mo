// Hybrid example: a decaying state whose gain switches at t = 1 and whose value
// is doubled at every sample instant t = 2, 3, ... (time events only).
model HybridTimeEvents
  discrete Real a(start = 1);
  Integer k(start = 0);
  Real x(start = 1);
equation
  der(x) = -a*x;
  when time >= 1 then
    a = 3;
  end when;
  when sample(2, 1) then
    k = pre(k) + 1;
    reinit(x, 2*x);
  end when;
end HybridTimeEvents;
