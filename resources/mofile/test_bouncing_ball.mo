// Bouncing ball: a state event (h < 0) reinitializes the velocity with a restitution coefficient.
// Bounces accumulate (Zeno behavior) at about t = 1.35 for e = 0.5, so simulate up to t = 1.2 or so.
model BouncingBall
  parameter Real g = 9.81;
  parameter Real e = 0.5;
  Real h(start = 1);
  Real v(start = 0);
equation
  der(h) = v;
  der(v) = -g;
  when h < 0 then
    reinit(v, -e*pre(v));
  end when;
end BouncingBall;
