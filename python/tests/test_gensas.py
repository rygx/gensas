"""Tests of the Python interface. Run with: python3 -m unittest discover -s python/tests
(needs the shared library, see python/gensas/__init__.py)."""

import math
import os
import sys
import tempfile
import threading
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gensas  # noqa: E402


def job(source, end_time):
    return {
        "gensas": {"format": "job", "version": 1},
        "kind": "modelica",
        "model": {"source": source},
        "solver": {"endTime": end_time},
    }


DECAY = "model D\n  parameter Real k = 1;\n  Real x(start = 1);\nequation\n  der(x) = -k*x;\nend D;\n"
GAIN = (
    "model M\n  discrete Real a(start = 1);\n  Real x(start = 1);\nequation\n  der(x) = -a*x;\n"
    "  when time >= 1 then\n    a = 3;\n  end when;\nend M;\n"
)
BALL = (
    "model Ball\n  Real h(start = 1);\n  Real v(start = 0);\nequation\n  der(h) = v;\n  der(v) = -9.81;\n"
    "  when h < 0 then\n    reinit(v, -0.5*pre(v));\n  end when;\nend Ball;\n"
)


class GenSASTest(unittest.TestCase):
    def test_check_and_diagnostics(self):
        with gensas.Simulation(job(DECAY, 1.0)) as sim:
            self.assertEqual(sim.check(), [])
            sim.run()
        bad = "model M\n  Real x(start = 1);\nequation\n  der(x) = -y;\nend M;\n"
        with gensas.Simulation(job(bad, 1.0)) as sim:
            found = sim.check()
            self.assertEqual(len(found), 1)
            self.assertEqual((found[0].severity, found[0].line, found[0].source), ("error", 4, "checker"))
            self.assertIn("undeclared identifier 'y'", found[0].message)
            self.assertEqual(found[0].as_dict()["line"], 4)
            with self.assertRaises(gensas.GenSASError) as ctx:
                sim.run()
            self.assertEqual(ctx.exception.diagnostics[0].line, 4)
        warn = "model M\n  Real x(start = 1, nominal = 2);\nequation\n  der(x) = -x;\nend M;\n"
        with gensas.Simulation(job(warn, 1.0)) as sim:
            found = sim.check()
            self.assertEqual([d.severity for d in found], ["warning"])

    def test_model_info(self):
        src = ("model Tank \"a tank\"\n  parameter Real k(unit = \"1/s\") = 2 \"rate\";\n"
               "  Real x(start = 1, unit = \"m\");\n  Real y;\nequation\n  der(x) = -k * x;\n  y = 2 * x;\nend Tank;\n")
        with gensas.Simulation(job(src, 1.0)) as sim:
            self.assertIsNone(sim.model_info())
            self.assertEqual(sim.check(), [])
            info = sim.model_info()
            self.assertEqual((info["name"], info["description"], info["hybrid"]), ("Tank", "a tank", False))
            kinds = {v["name"]: v["kind"] for v in info["variables"]}
            self.assertEqual(kinds, {"x": "state", "y": "algebraic"})
            x = [v for v in info["variables"] if v["name"] == "x"][0]
            self.assertEqual((x["start"], x["unit"]), (1.0, "m"))
            k = info["parameters"][0]
            self.assertEqual((k["name"], k["value"], k["unit"], k["constant"]), ("k", 2.0, "1/s", False))

    def test_model_selection(self):
        src = ("model A\n  Real a(start = 1);\nequation\n  der(a) = -a;\nend A;\n"
               "model B\n  Real b(start = 1);\nequation\n  der(b) = -b;\nend B;\n")
        with gensas.Simulation(job(src, 1.0)) as sim:
            self.assertEqual(sim.check(), [])
            info = sim.model_info()
            self.assertEqual((info["name"], info["models"], info["complete"]), ("B", ["A", "B"], True))
        chosen = job(src, 1.0)
        chosen["model"]["className"] = "A"
        self.assertIn("a", gensas.run(chosen).variables)
        chosen["model"]["className"] = "C"
        with gensas.Simulation(chosen) as sim:
            found = sim.check()
            self.assertIn("no model named 'C'", found[0].message)
            info = sim.model_info()
            self.assertEqual((info["models"], info["complete"]), (["A", "B"], False))

    def test_parameter_overrides(self):
        src = "model P\n  parameter Real k = 1;\n  Real x(start = 1);\nequation\n  der(x) = -k * x;\nend P;\n"
        job_json = job(src, 1.0)
        job_json["model"]["parameters"] = {"k": 3}
        job_json["model"]["start"] = {"x": 2}
        r = gensas.run(job_json)
        self.assertAlmostEqual(r.variables["x"][-1], 2 * math.exp(-3.0), places=4)
        job_json["model"]["parameters"] = {"nope": 1}
        with gensas.Simulation(job_json) as sim:
            found = sim.check()
            self.assertIn("parameters.nope", found[0].message)

    def test_versions(self):
        self.assertTrue(gensas.version())
        self.assertGreaterEqual(gensas.api_version(), 2)

    def test_decay(self):
        r = gensas.run(job(DECAY, 2.0))
        self.assertAlmostEqual(r.time[-1], 2.0, places=9)
        self.assertEqual(list(r.variables), ["x"])
        self.assertEqual(len(r.variables["x"]), len(r.time))
        self.assertAlmostEqual(r.variables["x"][-1], math.exp(-2.0), places=4)
        self.assertEqual(r.events, [])

    def test_hybrid_events(self):
        r = gensas.run(job(GAIN, 2.0))
        self.assertEqual(len(r.events), 1)
        self.assertAlmostEqual(r.events[0]["time"], 1.0, places=9)
        self.assertEqual(r.events[0]["assignments"], {"a": 3.0})
        self.assertAlmostEqual(r.discrete["a"][-1], 3.0)
        self.assertAlmostEqual(r.variables["x"][-1], math.exp(-1.0) * math.exp(-3.0), places=4)
        ball = gensas.run(job(BALL, 1.2))
        self.assertEqual(len(ball.events), 3)
        self.assertIn("v", ball.events[0]["reinits"])

    def test_terminate_and_assert(self):
        stop = ("model M\n  Real x(start = 1);\nequation\n  der(x) = -x;\n"
                "  when x < 0.3 then\n    terminate(\"done\");\n  end when;\nend M;\n")
        r = gensas.run(job(stop, 5.0))
        self.assertAlmostEqual(r.termination["time"], math.log(1.0 / 0.3), places=5)
        self.assertEqual(r.termination["message"], "done")
        self.assertAlmostEqual(r.time[-1], r.termination["time"], places=9)
        self.assertIsNone(gensas.run(job(DECAY, 1.0)).termination)
        failing = "model M\n  Real x(start = 1);\nequation\n  der(x) = -x;\n  assert(x > 0.5, \"too low\");\nend M;\n"
        with self.assertRaises(gensas.GenSASError) as ctx:
            gensas.run(job(failing, 2.0))
        self.assertIn("too low", str(ctx.exception))
        warning = failing.replace("\"too low\"", "\"too low\", AssertionLevel.warning")
        with gensas.Simulation(job(warning, 2.0)) as sim:
            sim.run()
            self.assertTrue(any("too low" in w for w in sim.warnings))

    def test_errors(self):
        with self.assertRaises(gensas.GenSASError) as ctx:
            gensas.Simulation({"gensas": {"format": "job", "version": 1}, "kind": "modelica", "model": {"source": "x"}, "colour": 1})
        self.assertIn("unknown key colour", str(ctx.exception))
        with self.assertRaises(gensas.GenSASError) as ctx:
            gensas.run(job("model M\n  Real x;\nequation\n  der(x) = \n", 1.0))
        self.assertIn("error:", str(ctx.exception))  # the parser's message is included

    def test_progress_and_cancel(self):
        seen = []
        gensas.run(job(DECAY, 5.0), progress=lambda t, end: seen.append((t, end)) and False)
        self.assertGreaterEqual(len(seen), 2)
        self.assertTrue(all(abs(end - 5.0) < 1e-12 for _, end in seen))
        self.assertEqual([t for t, _ in seen], sorted(t for t, _ in seen))
        with self.assertRaises(gensas.Cancelled):
            gensas.run(job(DECAY, 5.0), progress=lambda t, end: True)

        def boom(t, end):
            raise ValueError("from callback")

        with self.assertRaises(ValueError):
            gensas.run(job(DECAY, 5.0), progress=boom)

    def test_job_file_and_writers(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "m.mo"), "w") as f:
                f.write(DECAY)
            with open(os.path.join(d, "job.json"), "w") as f:
                f.write('{"gensas":{"format":"job","version":1},"kind":"modelica","model":{"file":"m.mo"},"solver":{"endTime":1}}')
            import pathlib

            with gensas.Simulation(pathlib.Path(d) / "job.json", base_dir=pathlib.Path(d)) as sim:
                r = sim.run()
                self.assertAlmostEqual(r.variables["x"][-1], math.exp(-1.0), places=4)
                for fmt in ("json", "csv", "mat"):
                    out = os.path.join(d, "out." + fmt)
                    sim.write(fmt, out)
                    self.assertGreater(os.path.getsize(out), 0)
                with self.assertRaises(gensas.GenSASError):
                    sim.write("xml", os.path.join(d, "out.xml"))

    def test_warnings(self):
        with gensas.Simulation(job(DECAY, 5000.0)) as sim:
            self.assertEqual(len(sim.warnings), 1)
            self.assertIn("limited", sim.warnings[0])
            self.assertAlmostEqual(sim.run().time[-1], 1000.0, places=6)

    def test_closed_simulation(self):
        sim = gensas.Simulation(job(DECAY, 1.0))
        sim.close()
        sim.close()  # closing twice is fine
        with self.assertRaises(gensas.GenSASError):
            sim.run()
        with self.assertRaises(gensas.GenSASError):
            sim.write("json", "unused.json")

    def test_threads(self):
        results = {}

        def work(k):
            src = DECAY.replace("k = 1", "k = %g" % (0.5 + 0.25 * k))
            results[k] = gensas.run(job(src, 2.0)).variables["x"][-1]

        threads = [threading.Thread(target=work, args=(k,)) for k in range(6)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        for k in range(6):
            self.assertAlmostEqual(results[k], math.exp(-2.0 * (0.5 + 0.25 * k)), places=4)


if __name__ == "__main__":
    unittest.main()
