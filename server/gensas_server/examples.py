"""Example models offered by GET /examples. Each one runs on its own (tests check that)."""

EXAMPLES = [
    {
        "id": "decay",
        "title": "Exponential decay",
        "description": "One state with a parameter: der(x) = -k*x.",
        "endTime": 5.0,
        "source": """model Decay "Exponential decay"
  parameter Real k = 1 "decay rate";
  Real x(start = 1) "amount";
equation
  der(x) = -k * x;
end Decay;
""",
    },
    {
        "id": "oscillator",
        "title": "Damped oscillator",
        "description": "Mass on a spring with friction; the velocity is a second state.",
        "endTime": 10.0,
        "source": """model Oscillator "Damped oscillator"
  parameter Real m = 1 "mass";
  parameter Real c = 0.2 "damping";
  parameter Real k = 4 "spring constant";
  Real x(start = 1, unit = "m") "position";
  Real v(start = 0, unit = "m/s") "velocity";
equation
  der(x) = v;
  m * der(v) = -c * v - k * x;
end Oscillator;
""",
    },
    {
        "id": "bouncing_ball",
        "title": "Bouncing ball",
        "description": "A state event (the ball reaches the floor) reverses the velocity. The ball loses half its speed at each bounce, so the bounces pile up at about 1.35 s (Zeno behavior); the end time stays before that.",
        "endTime": 1.2,
        "source": """model BouncingBall "Bouncing ball"
  parameter Real g = 9.81;
  parameter Real e = 0.5 "restitution";
  Real h(start = 1) "height";
  Real v(start = 0) "velocity";
equation
  der(h) = v;
  der(v) = -g;
  when h < 0 then
    reinit(v, -e * pre(v));
  end when;
end BouncingBall;
""",
    },
    {
        "id": "pendulum",
        "title": "Pendulum",
        "description": "Nonlinear: the angle enters through sin().",
        "endTime": 10.0,
        "source": """model Pendulum "Pendulum"
  parameter Real g = 9.81;
  parameter Real l = 1 "length";
  Real theta(start = 1) "angle";
  Real omega(start = 0) "angular velocity";
equation
  der(theta) = omega;
  der(omega) = -g / l * sin(theta);
end Pendulum;
""",
    },
    {
        "id": "predator_prey",
        "title": "Predator and prey",
        "description": "Lotka-Volterra: two coupled states with a product term.",
        "endTime": 15.0,
        "source": """model PredatorPrey "Lotka-Volterra"
  parameter Real a = 1.1;
  parameter Real b = 0.4;
  parameter Real c = 0.4;
  parameter Real d = 0.1;
  Real prey(start = 10);
  Real predator(start = 10);
equation
  der(prey) = a * prey - b * prey * predator;
  der(predator) = d * prey * predator - c * predator;
end PredatorPrey;
""",
    },
    {
        "id": "step_input",
        "title": "Step in the input",
        "description": "A time event: the gain changes at t = 1.",
        "endTime": 3.0,
        "source": """model StepInput "First order system with a step"
  Real x(start = 0);
  Real u;
equation
  u = if time < 1 then 0 else 1;
  der(x) = u - x;
end StepInput;
""",
    },
    {
        "id": "two_tanks",
        "title": "Two tanks (components)",
        "description": "A class used twice as a component; the second tank is fed by the first.",
        "endTime": 15.0,
        "source": """model Tank "A tank that drains through a valve"
  parameter Real k = 0.5 "outflow coefficient";
  Real h(start = 1) "level";
  input Real qin "inflow";
  Real qout "outflow";
equation
  der(h) = qin - qout;
  qout = k * h;
end Tank;

model TwoTanks "Two tanks in a row"
  Tank a(h(start = 2), k = 0.4), b(h(start = 0), k = 0.2);
equation
  a.qin = 0;
  b.qin = a.qout;
end TwoTanks;
""",
    },
    {
        "id": "rc_circuit",
        "title": "RC circuit (connectors)",
        "description": "A source, a resistor and a capacitor joined with connect(); the capacitor charges with time constant R*C.",
        "endTime": 5.0,
        "source": """connector Pin
  Real v;
  flow Real i;
end Pin;
model Resistor
  parameter Real R = 1;
  Pin p, n;
  Real u;
equation
  u = p.v - n.v;
  0 = p.i + n.i;
  p.i = u / R;
end Resistor;
model Capacitor
  parameter Real C = 1;
  Pin p, n;
  Real u(start = 0);
equation
  u = p.v - n.v;
  0 = p.i + n.i;
  C * der(u) = p.i;
end Capacitor;
model Source
  parameter Real V = 1;
  Pin p, n;
equation
  V = p.v - n.v;
  0 = p.i + n.i;
end Source;
model Ground
  Pin p;
equation
  p.v = 0;
end Ground;
model RC "RC circuit"
  Source s(V = 5);
  Resistor r(R = 2);
  Capacitor c(C = 0.5);
  Ground g;
equation
  connect(s.p, r.p);
  connect(r.n, c.p);
  connect(c.n, s.n);
  connect(s.n, g.p);
end RC;
""",
    },
    {
        "id": "chua",
        "title": "Chua circuit (chaos)",
        "description": "Three states and a piecewise linear diode (abs): the double-scroll attractor. A start value changed in the fourth digit gives a different curve after about 30 s. Plot x against z for the attractor.",
        "endTime": 100.0,
        "source": """model Chua "Chua circuit, dimensionless form"
  parameter Real alpha = 15.6;
  parameter Real beta = 28;
  parameter Real m0 = -1.143 "slope of the diode in the middle segment";
  parameter Real m1 = -0.714 "slope of the diode in the outer segments";
  Real x(start = 0.7) "voltage across the first capacitor";
  Real y(start = 0) "voltage across the second capacitor";
  Real z(start = 0) "inductor current";
  Real h "current through the diode";
equation
  h = m1 * x + 0.5 * (m0 - m1) * (abs(x + 1) - abs(x - 1));
  der(x) = alpha * (y - x - h);
  der(y) = x - y + z;
  der(z) = -beta * y;
end Chua;
""",
    },
]
