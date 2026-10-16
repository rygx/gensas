# Supported Modelica

ModelicaSAS simulates one flat Modelica `model`. This page lists what is accepted; anything else is rejected with an error that names the line, so an unsupported construct never gives a silently wrong result. For the command line see [Usage](Usage.md), for job files [Job format](JobFormat.md).

## Model structure

```modelica
model Name "description"
  parameter Real k = 2;
  Real x(start = 1), y;
initial equation
  x = 1;
equation
  der(x) = -k*x;
  y = sin(x);
end Name;
```

* A file can hold several `model`s, `function`s and `package`s. The last model is simulated; `model.className` in a job (see [Job format](JobFormat.md)) picks another. `//` and `/* */` comments, descriptive strings and `annotation(...)` clauses are accepted and ignored.
* Not supported: `extends`, nested models, `connector`/`connect` and `initial algorithm` sections (they give an error). See [Components](#components) for models used inside models.

## Declarations

* Types `Real`, `Integer` and `Boolean`, with the prefixes `parameter`, `constant`, `discrete`, `input` and `output`. Several variables can share a declaration (`Real x(start = 1), z;`).
* Types of your own: `type Volt = Real(unit = "V", start = 5);`, then `Volt u;` or `parameter Volt limit = 12;`. A type is made from `Real`, `Integer`, `Boolean` or a type defined before it (`type Level = Volt(start = 2);`), can stand in a package (`Phys.Rate k;`) and may stand before or after the model. Its modifications are applied to every variable of the type, and the modifications written at the variable come after and win. Array types (`type V3 = Real[3];`) are not supported; write the size at the variable (`Volt u[3];`).
* Modifications: `start` (a constant expression; it is the initial value of a state or the starting guess of an algebraic variable), `unit` and `fixed`. Other attributes are ignored with a warning.
* Bindings: `parameter Real p = 2*k;` (a constant expression, which may use parameters declared earlier) and `Real y = expr;` (the equation `y = expr`, which may use variables declared later). `input` needs a binding. Bindings on `discrete`, `Integer` and `Boolean` variables are not supported; assign them in a `when` equation.
* Arrays are sets of scalar variables: `Real x[3];`, `Real[2,3] A;`, `parameter Real k[3] = {1, 2, 3};`, `Real x[n](each start = 1)` or `start = {1, 2, 3}`. Subscripts are resolved when the model is read, so they must be constants: numbers, parameters declared earlier or loop variables (`x[i+1]`). Elements are named `x[2]` and `A[1,3]` in the results. `size(x, k)`, `sum(x)` and `prod(x)` are available.
* Whole-array equations and assignments are expanded element by element: `x = y + 2*z;`, `der(x) = -x;`, `w = if s > 1 then x else -x;`, `c = {{1, 2}, {3, 4}}*s;`, `y := sin(x);` in an algorithm, and `reinit(x, {5, 6});`. Scalars are broadcast inside an expression, `+`, `-`, a product or quotient with a scalar, unary minus and functions are elementwise, and both sides of an equation must have the same size. Not supported: the product of two arrays (a matrix product; use a `for` loop and `sum`), `^` and division by an array, arrays as `when`/`if` conditions, and arrays in declaration bindings (`Real y[3] = x;`; write an equation).

## Equations

* `der(x) = expr;` for states and `lhs = rhs;` for algebraic equations. Algebraic equations may be implicit (`y*y = x`); they are solved with Newton's method from the start values, so give the unknowns start values near the solution. A singular Jacobian is reported as an error.
* `initial equation` sections set start values with `x = <constant>;` and `der(x) = <constant>;`.
* `for i in 1:n loop ... end for;` (also `first:step:last`, and several iterators `for i in 1:n, j in 1:m loop`) is unrolled in `equation` and `initial equation` sections and inside `if` and `when` branches.
* The number of states must equal the number of differential equations, and the number of algebraic variables the number of algebraic equations; otherwise the model is rejected.

## Packages

```modelica
package Phys
  constant Real g = 9.81;
  function drag
    input Real v;
    output Real f;
  algorithm
    f := 0.1*v*abs(v);
  end drag;
end Phys;

model Fall
  Real v(start = 0);
equation
  der(v) = Phys.g - Phys.drag(v);
end Fall;
```

A `package` holds `constant`s and `function`s (and packages nested in it), used with qualified names such as `Phys.g` and `Phys.drag(v)`. A constant is a number or a constant expression of earlier constants of the same package. Functions of a package can call each other without the prefix. A package can also hold types (see Declarations). Not supported: models and `import` in a package.

## Components

```modelica
model Tank
  parameter Real k = 0.5;
  Real h(start = 1);
  input Real qin;
  Real qout;
equation
  der(h) = qin - qout;
  qout = k*h;
end Tank;

model TwoTanks
  Tank a(h(start = 2), k = 0.4), b(h(start = 0), k = 0.2);
equation
  a.qin = 0;
  b.qin = a.qout;
end TwoTanks;
```

A `model`, `block` or `record` of the source can be used as a component of the simulated model. Its variables, parameters and equations are added to the model once for each component, with the component's name in front: `a.h`, `b.h`, `a.k` (that is also how they are named in the results, and how a job names them in `model.parameters` and `model.start`). Components can contain components. In expressions `a.h` is the variable, and `a.h[2]` an element of an array in the component.

* Modifications: `Tank a(k = 0.4)` sets a parameter (a constant expression, which may use parameters declared earlier), `a(h(start = 2))` sets a start value, `a(inner(k = 1))` or `a(inner.k = 1)` reach into a component of the component. A modification of a variable that is not a parameter (`Tank a(qout = 1)`) adds the equation `a.qout = 1`. The job's values win over the modifications, and a name the class does not have is an error.
* `parameter Pt c(x = 2, y = 3);` makes every variable of the record `Pt` a parameter, so records are the way to group parameters. `input` variables of a component can stay without a binding: the model around it must give them an equation (`b.qin = a.qout;`).
* Names of packages, functions, types and classes, and `time`, are not renamed inside a component. A class that contains itself is an error.
* Not supported: arrays of components (`Tank t[3];`), functions defined inside a class that is used as a component (define them before the class), `connector`s and `connect` equations, `extends`, redeclarations and `inner`/`outer`.

## Algorithm sections

```modelica
algorithm
  y := 2*x;
  if x < 1 then
    z := log(x);
  else
    z := 0;
  end if;
  for i in 1:n loop
    w := w + v[i];   // after w := 0;
  end for;
```

An `algorithm` section is turned into equations, one for each assigned variable, so everything that works in an `equation` section works on its results. Statements: assignments `v := expr;`, `if`/`elseif`/`else`, `for` loops (unrolled, as in equations), `when` statements, `assert(...)` and `reinit(...)` in a `when` body.

* A variable can be assigned several times; later statements see the latest value (`s := y + 1; s := s*0.5;`).
* An `if` statement becomes an if-expression, so the branch that is not taken is not evaluated (`log(x)` above is only evaluated where `x < 1`), and a condition on continuous variables is a state event. A variable assigned in only some branches must have a value before the `if`.
* A variable must not be read before it is assigned in the section: its previous value is not available (`y := y + 1;` is an error). In a `when` body this is allowed for discrete variables: `cnt := cnt + 1;` means `cnt := pre(cnt) + 1;`.
* Not supported: `while` loops, `break` and `return`, `when` and `assert` inside an `if` statement, `if` inside a `when` body, `initial algorithm`.

## Functions

```modelica
function poly
  input Real a;
  input Real b = 2;
  output Real y;
protected
  Real t = a*a;
algorithm
  if t > 3 then
    y := t - b;
  else
    y := t + b;
  end if;
end poly;
```

A function can be defined before the `model` or inside it, and is called as `poly(x)` or `poly(x, 5)`. The body is an `algorithm` section (see above) and the call is replaced by the result expression, so a function costs nothing at run time and its `if` branches are state events like any other.

* One output, scalar `Real`, `Integer` and `Boolean` inputs and protected variables, default values for inputs.
* A function can call other functions. It can use parameters and constants of the model, but no other model variable.
* A call can be used in equations, `when` and `if` conditions, parameter values and start values (a function used in a parameter value must be defined before the parameter).
* Not supported: recursion, several outputs, array inputs and outputs, named arguments (`f(a = 1)`), `external` functions, `when` and `assert` in the body.

## Expressions

* Operators `+ - * / ^`, unary minus, relations `< <= > >= == <>` and `and`, `or`, `not`.
* A small integer power such as `x^2` or `x^-1` becomes products, `x^0.5` becomes `sqrt`, and any other power becomes `exp(y*log(x))`, which needs `x > 0`.
* Functions: `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `exp`, `log`, `log10`, `sqrt`, and `abs`, `min`, `max`, `sign`. A function of constants (`exp(1)`) is evaluated when the model is read.
* `abs`, `min`, `max` and `sign` are turned into if-expressions, and their switching points are state events, so the solution is exact. A variable that hovers around a switching point (`sign(v)` friction at `v = 0`) makes the events pile up and stops the run with an error; `--piecewise smooth` replaces the functions by smooth approximations instead (see [Usage](Usage.md)).
* `atan2` is rewritten into `atan` with if-branches; the switches at `x = 0` and at the branch cut are state events.
* `noEvent(e)` is accepted and treated like `e`.
* `if c then a else b` expressions, `time`, and `pre(d)` for discrete variables.

## Hybrid models

* Discrete variables (`discrete Real`, `Integer`, `Boolean`) change only at events.
* `if` expressions and `if` equations, with conditions on parameters, discrete variables, `time` or continuous variables. An `if` equation may change which equations are active (a variable can be a state in one mode and algebraic in another) as long as every mode involves the same set of variables.
* `when` equations with `elsewhen`, with conditions
  * `time >= e` and `sample(start, interval)` (time events),
  * relations on continuous variables, such as `h < 0` (state events, located along the solution by bisection),
  * `initial()`, which runs once at the start time,
  * `edge(b)` and `change(d)` on a single discrete variable.

  The body assigns discrete variables (`d = expr;`) and re-initializes states with `reinit(x, expr);`.
* `terminate("message");` in a `when` equation ends the run successfully at that event; the result stops there and records the message.
* `assert(condition, "message");` in the equation section is checked along the solution (a state event when the condition turns false) and in `when` bodies at the event. A violation stops the run with an error naming the time; with `AssertionLevel.warning` as a third argument it is reported once as a warning and the run goes on. `terminate` and `assert` are not supported inside `if` equations.
* A simulation whose events accumulate in finite time (Zeno behavior, for example a bouncing ball) stops with an error.

## What the solver checks

The solver is a semi-analytical series method, so a wrong result would be easy to overlook; these checks turn the common causes into errors or warnings.

* **Function domains.** The arguments of `log`, `sqrt`, `asin`, `acos` and the denominators of divisions are checked along every step. A model that leaves a domain, divides by zero or passes a pole (a denominator that changes sign) stops with an error that names the time and the operand. A guard such as `if x > 0 then log(x) else 0` keeps the argument inside the domain; in hybrid models the check applies to the part of the solution that is kept after an event. A pole where the denominator does not change sign, such as `1/u^2`, is not detected.
* **No progress.** If the solver cannot advance within the error tolerance, the run fails with the time of the stall instead of returning a shortened result.
* **End time.** A run that finishes always ends exactly at the requested end time.
* **Warnings.** A step that starts with a residual above the error tolerance is reported as a warning in the result (at most 10 per run), and so are failed `assert` calls with `AssertionLevel.warning` and out-of-range solver options that were clamped.

## Examples

`resources/mofile/` has runnable models: `test_solve_ode.mo` (an ODE system), `test_hybrid_time_events.mo` (time events), `test_bouncing_ball.mo` (a state event with `reinit`) and `test_array_cascade.mo` (arrays and a `for` loop).
