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
        "description": "A state event (the ball reaches the floor) reverses the velocity.",
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
]
