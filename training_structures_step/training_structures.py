# -*- coding: utf-8 -*-

"""Non-graphical part of the Training Structures step in a SEAMM flowchart"""

import logging
from pathlib import Path
import importlib.resources
import pprint  # noqa: F401

import numpy as np

import training_structures_step
from training_structures_step import generate
import molsystem
import seamm
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Training Structures")

# Add this module's properties to the standard properties
path = importlib.resources.files("training_structures_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


class TrainingStructures(seamm.Node):
    """
    The non-graphical part of a Training Structures step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : TrainingStructuresParameters
        The control parameters for Training Structures.

    See Also
    --------
    TkTrainingStructures,
    TrainingStructures, TrainingStructuresParameters
    """

    def __init__(
        self, flowchart=None, title="Training Structures", extension=None, logger=logger
    ):
        """A step for Training Structures in a SEAMM flowchart.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Training Structures {self}")

        super().__init__(
            flowchart=flowchart,
            title="Training Structures",
            extension=extension,
            module=__name__,
            logger=logger,
        )  # yapf: disable

        self._metadata = training_structures_step.metadata
        self.parameters = training_structures_step.TrainingStructuresParameters()

    @property
    def version(self):
        """The semantic version of this module."""
        return training_structures_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return training_structures_step.__git_revision__

    # ------------------------------------------------------------------ #
    # Description
    # ------------------------------------------------------------------ #
    def description_text(self, P=None):
        """Create the text description of what this step will do."""
        if not P:
            P = self.parameters.values_to_dict()

        n = P["number of samples"]
        if str(n).strip().lower() == "in proportion to the modes":
            count = (
                f"{P['samples per mode']} normal-mode samples per mode above "
                f"{P['lowest frequency']} (at least {P['minimum number of samples']}, "
                f"at most {P['maximum number of samples']})"
            )
        else:
            count = f"{n} normal-mode samples"
        text = (
            "Distorted structures of the current molecule, for training a "
            "machine-learned potential, from the Hessian of the model chemistry "
            f"defined earlier in the flowchart: {count}, each scaled to a harmonic "
            f"energy between {P['minimum target energy']} and "
            f"{P['maximum target energy']}"
        )
        if str(P["number of rotors"]) != "0":
            text += (
                f"; relaxed torsion scans in steps of {P['scan step']} of up to "
                f"{P['number of rotors']} heavy rotors at {P['generator level']}"
            )
        if P["bond stretches"] == "yes":
            text += (
                "; and two stretches of the median bond of each bonded element pair "
                f"({P['bond pairs']}), to a harmonic energy of {P['stretch energy']} "
                f"but at most {P['maximum stretch']}"
            )
        text += (
            ". A diatomic gets a stretch grid of "
            f"{P['diatomic grid points']} points instead."
        )
        checks = []
        if P["connectivity check"] == "yes":
            checks.append("keeps the minimum's bonds")
        if not self._is_none(P["prescreen energy"]):
            checks.append(
                f"lies within {P['prescreen energy']} of the minimum at "
                f"{P['generator level']}"
            )
        if checks:
            text += " Every point must pass a check that it " + " and ".join(checks)
            text += "."
        route = P.get("hessian route", "warm engine over MDI")
        if route == "separate tasks":
            text += (
                " The Hessian is central finite differences of the gradient, run as "
                "separate tasks."
            )
        else:
            text += f" The Hessian comes from the Model Chemistry ({route})."
            if P["check the hessian"] != "never":
                text += (
                    " It is checked against finite differences of the gradient "
                    f"({P['check the hessian']}) and replaced by a finite-difference "
                    "Hessian if it fails."
                )
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    # ------------------------------------------------------------------ #
    # Run
    # ------------------------------------------------------------------ #
    def run(self):
        """Make the training structures of the current molecule."""

        next_node = super().run(printer)
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        printer.important(__(self.description_text(P), indent=self.indent))
        printer.important("")

        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        system_db = self.get_variable("_system_db")
        system, reference = self.get_system_configuration(None)
        if reference.periodicity != 0:
            raise ValueError("Training Structures handles molecules, not crystals.")

        seed = self._seed_of(P["random seed"])
        rng = np.random.default_rng(seed)

        symbols = list(reference.atoms.symbols)
        n_atoms = len(symbols)
        if n_atoms < 2:
            raise ValueError("Training Structures needs at least two atoms.")
        x0 = np.asarray(reference.atoms.get_coordinates(fractionals=False), float)
        masses = np.asarray(reference.atoms.atomic_masses, dtype=float)
        factor = float(P["connectivity factor"])
        ref_bonds = generate.connectivity(symbols, x0, factor)
        table_bonds, orders, conjugated = self._bond_table(reference)
        if table_bonds and set(table_bonds) != ref_bonds:
            printer.important(
                __(
                    "Note: the bonds perceived at the minimum differ from the "
                    "structure's bond table; points are compared with the perceived "
                    "bonds.",
                    indent=4 * " ",
                )
            )

        driver = self._driver(P)
        if driver.model_chemistry() is None:
            raise ValueError(
                "Training Structures needs a Model Chemistry step before it."
            )

        data = {"random seed": seed}
        # 1. The Hessian at the Model Chemistry, checked
        H, data["hessian"] = self._hessian(P, driver, reference, x0, masses, rng)
        # What the engine actually ran for the Hessian (its calculation record, from
        # the program's own input and output, which stay in the engine's scratch)
        settings = engine_settings(driver.fingerprints)
        if settings:
            data["hessian"]["settings"] = settings

        # 2. The points, with the generator level for the scans and pre-screen
        points = [{"type": "minimum", "x": x0, "name": "minimum"}]
        dropped = {}
        level = P["generator level"]
        prescreen = (
            None
            if self._is_none(P["prescreen energy"])
            else self._as_float(P["prescreen energy"], "kJ/mol")
        )
        rotors = []
        if n_atoms > 2 and int(P["number of rotors"]) > 0:
            rotors = generate.find_rotors(
                symbols,
                table_bonds or sorted(ref_bonds),
                orders,
                conjugated,
                max_rotors=int(P["number of rotors"]),
            )
        need_cheap = prescreen is not None or rotors
        cheap_context = (
            driver.cheap_session(reference, level) if need_cheap else _no_session()
        )
        with cheap_context as cheap:
            if need_cheap and cheap is None:
                raise RuntimeError(
                    f"The generator level {level} cannot run as a warm engine here."
                )
            e_ref = cheap.energy(x0, units="kJ/mol") if cheap is not None else None

            def check(x, energy=None, stretched=None):
                """(reason or None, generator-level energy above the minimum). The
                bond a stretch point lengthens on purpose (stretched, i < j) is
                not judged by the connectivity check; every other bond is."""
                if P["connectivity check"] == "yes":
                    change = generate.connectivity_change(
                        ref_bonds, symbols, x, factor, exempt=stretched
                    )
                    if change is not None:
                        return f"bonds changed: {change}", None
                de = None
                if cheap is not None:
                    e = cheap.energy(x, units="kJ/mol") if energy is None else energy
                    de = e - e_ref
                    if prescreen is not None and de > prescreen:
                        return (
                            f"{de:.1f} kJ/mol above the minimum at {level}",
                            de,
                        )
                return None, de

            if n_atoms == 2:
                points += self._diatomic(P, H, x0, check, dropped)
            else:
                points += self._samples(P, H, x0, masses, symbols, check, rng, dropped)
                for k, rotor in enumerate(rotors, start=1):
                    points += self._scan(
                        P, cheap, rotor, k, x0, symbols, check, dropped
                    )
                if P["bond stretches"] == "yes":
                    points += self._stretches(P, H, x0, symbols, check, dropped)

        # 3. The new system, one configuration per point
        name = self._system_name(P["system name"], system.name)
        out_sys = system_db.create_combined_system([reference], name=name)
        base = out_sys.configuration
        for k, point in enumerate(points):
            conf = base if k == 0 else out_sys.copy_configuration(base)
            conf.atoms.set_coordinates(point["x"], fractionals=False)
            conf.name = point["name"]
            self._put_property(conf, "point type", point["type"], "", "str")
            for key, (prop, units, _type) in _POINT_PROPERTIES.items():
                if point.get(key) is not None:
                    self._put_property(conf, prop, point[key], units, _type)
        self._put_property(base, "hessian", data["hessian"], "", "json")
        self._put_property(base, "random seed", seed, "", "int")
        system_db.system = out_sys
        out_sys.configuration = base.id

        counts = {}
        for point in points:
            counts[point["type"]] = counts.get(point["type"], 0) + 1
        data["number of points"] = len(points)
        for key in ("sample", "scan", "stretch", "grid"):
            data[f"number of {key} points"] = counts.get(key, 0)
        data["number of rotors"] = len(rotors)
        data["dropped points"] = dropped
        data["number of dropped points"] = sum(len(v) for v in dropped.values())
        data["system name"] = name

        self.analyze(data=data)
        self.store_results(configuration=base, data=data)
        return next_node

    def _driver(self, P):
        """seamm_exec's HessianDriver for the Model Chemistry, on the chosen route."""
        from seamm_exec.hessian import HessianDriver

        return HessianDriver(
            self,
            route=P["hessian route"],
            who="Training Structures",
            print_to=printer,
        )

    # ------------------------------------------------------------------ #
    # The Hessian
    # ------------------------------------------------------------------ #
    def _hessian(self, P, driver, configuration, x0, masses, rng):
        """The Hessian (E_h/bohr²) at the Model Chemistry and how it was made,
        checked against finite differences of the gradient when asked.

        Over a warm engine the check's gradients each start from the reference
        structure's orbitals, with a tight SCF (seamm_exec's
        ``EngineSession.check_hessian``). A refused Hessian is replaced by finite
        differences of the gradient, made the same way and with every point
        checked for a state change, and the replacement is checked again: if it
        fails too, the molecule is refused (science 2026-11-26)."""
        how = P["check the hessian"]
        open_shell = int(configuration.spin_multiplicity or 1) > 1
        draw = rng.uniform()  # always drawn, so the points do not depend on how
        check = how == "always" or (
            how == "open shells"
            and (open_shell or draw < float(P["closed-shell check fraction"]))
        )
        info = {"open shell": open_shell, "checked": check}
        limits = {
            "rms": self._as_float(P["check displacement"], "Å"),
            "tolerance": float(P["hessian tolerance"]),
            "max_frequency": self._as_float(P["maximum frequency"], "1/cm"),
        }
        with driver.session(configuration) as session:
            if session is None:
                _, H, _ = driver.energy_and_hessian(configuration)
                info["route"] = "tasks"
                fd = getattr(driver, "fd_info", None)
                if fd:
                    info["finite differences"] = fd
                if check:
                    info["check"] = "not possible: the model chemistry runs as tasks"
                return np.asarray(H, dtype=float), info
            H, route = session.hessian()
            info["route"] = route
            fd = getattr(session, "fd_info", None)
            if fd:
                info["finite differences"] = fd
            if not check:
                return np.asarray(H, dtype=float), info
            result = self._check_hessian(session, H, x0, masses, rng, limits)
            info["errors"] = [round(e, 6) for e in result["errors"]]
            info["max frequency"] = round(result["max frequency"], 1)
            if result.get("reference guess") is not None:
                info["check from reference guess"] = result["reference guess"]
            if result["refused"]:
                printer.important(
                    __(
                        f"The {route} Hessian was refused: {result['reason']}. It is "
                        "replaced by finite differences of the gradient.",
                        indent=4 * " ",
                    )
                )
                info["refused"] = result["reason"]
                H, route = session.hessian(finite_differences=True)
                info["route"] = route
                fd = getattr(session, "fd_info", None)
                if fd:
                    info["finite differences"] = fd
                # The replacement is checked too, and must pass
                again = self._check_hessian(session, H, x0, masses, rng, limits)
                info["replacement check"] = {
                    "errors": [round(e, 6) for e in again["errors"]],
                    "max frequency": round(again["max frequency"], 1),
                    "reason": again["reason"],
                }
                if again["refused"]:
                    raise RuntimeError(
                        f"The {route} Hessian that replaced it was refused too: "
                        f"{again['reason']}. The molecule is not used."
                    )
        return np.asarray(H, dtype=float), info

    @staticmethod
    def _check_hessian(session, H, x0, masses, rng, limits):
        """The gradient-difference check: the session's, from the reference
        orbitals with a tight SCF (seamm-exec >= 2026.11.26), else with its
        plain gradients."""
        check = getattr(session, "check_hessian", None)
        if check is not None:
            return check(H, rng, coords_ang=x0, **limits)
        return generate.check_hessian(H, session.gradient, x0, masses, rng, **limits)

    # ------------------------------------------------------------------ #
    # The kinds of points
    # ------------------------------------------------------------------ #
    def _samples(self, P, H, x0, masses, symbols, check, rng, dropped):
        """Energy-targeted normal-mode samples."""
        from normal_mode_sampling_step.normal_mode_sampling import (
            ZERO_MODE_THRESHOLD_CM,
            geometry_guard,
            normal_mode_analysis,
            sample_normal_modes,
        )

        freqs, modes = normal_mode_analysis(H, masses, x0)
        keep = freqs > ZERO_MODE_THRESHOLD_CM
        if not self._is_none(P["lowest frequency"]):
            keep &= freqs >= self._as_float(P["lowest frequency"], "1/cm")
        freqs, modes = freqs[keep], modes[:, keep]
        n = self._number_of_samples(P, freqs.size)
        if n == 0 or freqs.size == 0:
            return []
        limits = {
            key: (None if self._is_none(P[key]) else self._as_float(P[key], units))
            for key, units in (
                ("maximum atom displacement", "Å"),
                ("maximum bond change", "Å"),
                ("maximum angle change", "degree"),
            )
        }
        geometric = geometry_guard(
            x0,
            sorted(generate.connectivity(symbols, x0, float(P["connectivity factor"]))),
            max_displacement=limits["maximum atom displacement"],
            max_bond_change=limits["maximum bond change"],
            max_angle_change=limits["maximum angle change"],
        )
        energies = {}
        reasons = []

        def guard(dx):
            if geometric is not None:
                reason = geometric(dx)
                if reason is not None:
                    reasons.append(reason)
                    return reason
            reason, de = check(x0 + dx)
            if reason is not None:
                reasons.append(reason)
                return reason
            energies[dx.tobytes()] = de
            return None

        cap = None if self._is_none(P["amplitude cap"]) else float(P["amplitude cap"])
        displacements, targets, _, _, n_dropped = sample_normal_modes(
            freqs,
            modes,
            masses,
            n,
            distribution="energy",
            energy_range=(
                self._as_float(P["minimum target energy"], "kJ/mol"),
                self._as_float(P["maximum target energy"], "kJ/mol"),
            ),
            amplitude_cap=cap,
            guard=guard,
            max_draws=int(P["maximum draws"]),
            rng=rng,
        )
        if n_dropped:
            dropped["sample"] = [
                f"{n_dropped} of {n} samples failed all {P['maximum draws']} draws"
            ]
        self._redraws = len(reasons)
        return [
            {
                "type": "sample",
                "name": f"sample {k}",
                "x": x0 + dx,
                "target energy": float(e),
                "harmonic energy": float(e),
                "prescreen energy": energies.get(dx.tobytes()),
            }
            for k, (dx, e) in enumerate(zip(displacements, targets), start=1)
        ]

    def _scan(self, P, cheap, rotor, k, x0, symbols, check, dropped):
        """A relaxed scan of one rotor at the generator level, from the minimum
        around the circle, each point relaxed with the torsion fixed."""
        from ase import Atoms
        from ase.constraints import FixInternals
        from ase.data import atomic_numbers
        from ase.optimize import BFGS
        from seamm_exec.hessian import _MDICalculator

        a, i, j, b = rotor["dihedral"]
        label = f"{symbols[i]}{i + 1}-{symbols[j]}{j + 1}"
        start = generate.dihedral(x0, a, i, j, b)
        fmax = self._as_float(P["relaxation convergence"], "eV/Å")
        steps = int(P["relaxation max steps"])
        step = self._as_float(P["scan step"], "degree")
        numbers = [atomic_numbers[s] for s in symbols]
        points = []
        x = x0
        for target in generate.scan_angles(start, step):
            y = generate.set_dihedral(x, rotor, target)
            atoms = Atoms(numbers=numbers, positions=y)
            atoms.set_constraint(FixInternals(dihedrals_deg=[[target, [a, i, j, b]]]))
            atoms.calc = _MDICalculator(cheap.engine)
            converged = BFGS(atoms, logfile=None).run(fmax=fmax, steps=steps)
            x = np.asarray(atoms.positions, dtype=float)
            energy = atoms.get_potential_energy() * _EV_TO_KJ
            angle = target - start
            name = f"scan {k} {label} {angle:+.0f}"
            reason, de = check(x, energy=energy)
            if not converged:
                reason = f"the relaxation did not converge in {steps} steps"
            if reason is not None:
                dropped.setdefault("scan", []).append(f"{name}: {reason}")
                continue
            points.append(
                {
                    "type": "scan",
                    "name": name,
                    "x": x,
                    "scan rotor": [a + 1, i + 1, j + 1, b + 1],
                    "scan angle": float(angle),
                    "prescreen energy": de,
                }
            )
        return points

    def _stretches(self, P, H, x0, symbols, check, dropped):
        """Two stretches of the median bond of each bonded element pair."""
        pairs = None
        if str(P["bond pairs"]).strip().lower() != "all":
            pairs = {p for p in str(P["bond pairs"]).replace(";", ",").split(",") if p}
        bonds = sorted(
            generate.connectivity(symbols, x0, float(P["connectivity factor"]))
        )
        energy = self._as_float(P["stretch energy"], "kJ/mol") / generate.KJ_HARTREE
        maximum = self._as_float(P["maximum stretch"], "Å")
        points = []
        for cls, (i, j) in generate.bond_classes(symbols, x0, bonds, pairs).items():
            a = generate.stretch_amplitude(H, x0, i, j, energy, maximum)
            v = generate.stretch_vector(x0, i, j) / generate.BOHR
            k_eff = float(v @ H @ v)
            for sign in (1.0, -1.0):
                x = generate.stretch(x0, i, j, sign * a)
                name = f"stretch {cls} {'+' if sign > 0 else '-'}"
                reason, de = check(x, stretched=(i, j))
                if reason is not None:
                    dropped.setdefault("stretch", []).append(f"{name}: {reason}")
                    continue
                points.append(
                    {
                        "type": "stretch",
                        "name": name,
                        "x": x,
                        "bond class": cls,
                        "bond": [i + 1, j + 1],
                        "bond change": float(sign * a),
                        "harmonic energy": 0.5 * k_eff * a * a * generate.KJ_HARTREE,
                        "prescreen energy": de,
                    }
                )
        return points

    def _diatomic(self, P, H, x0, check, dropped):
        """A diatomic's stretch grid, in place of samples and stretches."""
        energy = (
            self._as_float(P["diatomic grid energy"], "kJ/mol") / generate.KJ_HARTREE
        )
        v = generate.stretch_vector(x0, 0, 1) / generate.BOHR
        k_eff = float(v @ H @ v)
        points = []
        grid = generate.diatomic_grid(H, x0, energy, int(P["diatomic grid points"]))
        for n, a in enumerate(grid, start=1):
            x = generate.stretch(x0, 0, 1, a)
            name = f"grid {n}"
            reason, de = check(x, stretched=(0, 1))
            if reason is not None:
                dropped.setdefault("grid", []).append(f"{name}: {reason}")
                continue
            points.append(
                {
                    "type": "grid",
                    "name": name,
                    "x": x,
                    "bond change": float(a),
                    "harmonic energy": 0.5 * k_eff * a * a * generate.KJ_HARTREE,
                    "prescreen energy": de,
                }
            )
        return points

    # ------------------------------------------------------------------ #
    # Output
    # ------------------------------------------------------------------ #
    def analyze(self, data=None, indent="", **kwargs):
        """Summarize the points made."""
        if data is None:
            return
        parts = [
            f"{data['number of sample points']} samples",
            f"{data['number of scan points']} scan points over "
            f"{data['number of rotors']} rotors",
            f"{data['number of stretch points']} bond stretches",
        ]
        if data["number of grid points"]:
            parts = [f"{data['number of grid points']} stretch-grid points"]
        text = (
            f"Made {data['number of points']} structures in the system "
            f"'{data['system name']}': the minimum, " + ", ".join(parts) + "."
        )
        printer.important(__(text, indent=4 * " "))
        hessian = data["hessian"]
        text = f"The Hessian is {hessian['route']}"
        if hessian.get("refused"):
            text += f" (the first one was refused: {hessian['refused']})"
        elif "errors" in hessian:
            text += (
                f"; checked against the gradient, largest error "
                f"{max(hessian['errors']):.4f}"
            )
        printer.important(__(text + ".", indent=4 * " "))
        if data["number of dropped points"]:
            printer.important(
                __(
                    f"{data['number of dropped points']} points were dropped:",
                    indent=4 * " ",
                )
            )
            for kind, items in data["dropped points"].items():
                # group the points refused for the same reason
                by_reason = {}
                for item in items:
                    name, _, reason = item.partition(": ")
                    by_reason.setdefault(reason or name, []).append(name)
                for reason, names in by_reason.items():
                    if len(names) == 1 or reason == names[0]:
                        text = f"{kind}: {names[0]}" + (
                            f": {reason}" if reason != names[0] else ""
                        )
                    else:
                        text = f"{kind}: {len(names)} points ({names[0]} ... "
                        text += f"{names[-1]}): {reason}"
                    printer.important(__(text, indent=8 * " ", wrap=False))
        printer.important(
            __(f"The random seed was {data['random seed']}.", indent=4 * " ")
        )
        printer.important("")

    def _put_property(self, conf, name, value, units, _type="float"):
        """Store a per-point property (the names in data/properties.csv)."""
        name = f"{name}#TrainingStructures#generation"
        properties = conf.properties
        if not properties.exists(name):
            properties.add(name, _type, units=units or None, noerror=True)
        properties.put(name, value)

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    @staticmethod
    def _bond_table(configuration):
        """(bonds, orders, conjugated) of a configuration's bond table, 0-based;
        conjugation from RDKit. Empty when it has no bonds."""
        if configuration.bonds.n_bonds == 0:
            return [], {}, set()
        ids = list(configuration.atoms.ids)
        index = {aid: k for k, aid in enumerate(ids)}
        table = configuration.bonds.get_as_dict()
        bonds, orders = [], {}
        for i, j, o in zip(table["i"], table["j"], table["bondorder"]):
            key = tuple(sorted((index[i], index[j])))
            bonds.append(key)
            orders[key] = int(o)
        conjugated = set()
        try:
            mol = configuration.to_RDKMol()
            for b in mol.GetBonds():
                if b.GetIsConjugated():
                    conjugated.add(
                        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
                    )
        except Exception as e:  # pragma: no cover
            logger.warning(f"No conjugation from RDKit: {e}")
        return bonds, orders, conjugated

    @staticmethod
    def _system_name(how, structure):
        """The new system's name: '<structure> training', '<structure> redraw', or
        the text given with {structure} replaced."""
        how = str(how)
        if how == "from structure":
            return f"{structure} training"
        if how == "for a redraw":
            # a redraw starts from the first pass's minimum, in '<x> training'
            if structure.endswith(" training"):
                structure = structure[: -len(" training")]
            return f"{structure} redraw"
        return how.replace("{structure}", structure)

    @staticmethod
    def _seed_of(seed):
        if isinstance(seed, str):
            if seed.strip() == "" or seed.strip().lower() == "random":
                return int(np.random.SeedSequence().entropy % 2**63)
            seed = int(seed)
        return int(seed)

    @staticmethod
    def _number_of_samples(P, n_modes):
        n = P["number of samples"]
        if str(n).strip().lower() != "in proportion to the modes":
            return int(n)
        wanted = int(np.floor(float(P["samples per mode"]) * n_modes + 0.5))
        low = int(P["minimum number of samples"])
        high = int(P["maximum number of samples"])
        return max(low, min(high, wanted))

    @staticmethod
    def _is_none(value):
        return isinstance(value, str) and value.strip().lower() in ("none", "")

    @staticmethod
    def _as_float(value, units):
        try:
            return value.m_as(units)
        except AttributeError:
            return float(value)


_EV_TO_KJ = 96.48533212331  # kJ/mol per eV

#: per-point values -> (property, units, type)
_POINT_PROPERTIES = {
    "target energy": ("target energy", "kJ/mol", "float"),
    "harmonic energy": ("harmonic energy", "kJ/mol", "float"),
    "prescreen energy": ("prescreen energy", "kJ/mol", "float"),
    "scan rotor": ("scan rotor", "", "json"),
    "scan angle": ("scan angle", "degree", "float"),
    "bond class": ("bond class", "", "str"),
    "bond": ("bond", "", "json"),
    "bond change": ("bond change", "Å", "float"),
}


def engine_settings(fingerprints):
    """The settings of the engine's Hessian run (else its last energy run), from
    the calculation records a warm engine reports ({"energy": record, "hessian":
    record}, see seamm_util.fingerprint.calculation_record): {"run", "fingerprint",
    "version", "keywords"}, the keywords being the input's '!' lines. None without
    records (e.g. the Hessian ran as separate tasks)."""
    if not fingerprints:
        return None
    for run in ("hessian", "energy"):
        record = fingerprints.get(run)
        if record:
            files = record.get("files") or {}
            text = "\n".join(str(v) for v in files.values())
            return {
                "run": run,
                "fingerprint": record.get("fingerprint"),
                "version": record.get("version"),
                "keywords": [
                    line.strip()
                    for line in text.splitlines()
                    if line.strip().startswith("!")
                ],
            }
    return None


class _no_session:
    """A context manager that gives no engine (no generator level needed)."""

    def __enter__(self):
        return None

    def __exit__(self, *args):
        return False
