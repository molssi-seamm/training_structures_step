# -*- coding: utf-8 -*-

"""Non-graphical part of the Export Training Set step in a SEAMM flowchart"""

import csv
import fnmatch
import logging
from pathlib import Path
import importlib.resources
import pprint  # noqa: F401
import re

import training_structures_step
from training_structures_step import export, generate
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
printer = printing.getPrinter("Export Training Set")

# Add this module's properties to the standard properties
path = importlib.resources.files("training_structures_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


class ExportTrainingSet(seamm.Node):
    """
    The non-graphical part of a Export Training Set step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : ExportTrainingSetParameters
        The control parameters for Export Training Set.

    See Also
    --------
    TkExportTrainingSet,
    ExportTrainingSet, ExportTrainingSetParameters
    """

    def __init__(
        self, flowchart=None, title="Export Training Set", extension=None, logger=logger
    ):
        """A step for Export Training Set in a SEAMM flowchart.

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
        logger.debug(f"Creating Export Training Set {self}")

        super().__init__(
            flowchart=flowchart,
            title="Export Training Set",
            extension=extension,
            module=__name__,
            logger=logger,
        )  # yapf: disable

        self._metadata = training_structures_step.export_metadata
        self.parameters = training_structures_step.ExportTrainingSetParameters()

    @property
    def version(self):
        """The semantic version of this module."""
        return training_structures_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return training_structures_step.__git_revision__

    def description_text(self, P=None):
        """Create the text description of what this step will do."""
        if not P:
            P = self.parameters.values_to_dict()
        text = (
            f"Export the labelled training points of {P['systems']}, one "
            "extended-XYZ file per molecule with energies and forces "
            f"({P['model chemistry']}) in {P['energy units']}, each energy also "
            "relative to the molecule's minimum."
        )
        if not self._is_none(P["energy window"]):
            text += f" Points more than {P['energy window']} above it are dropped."
        if not self._is_none(P["spin tolerance"]):
            text += (
                " For open shells, points whose <S**2> differs from the minimum's "
                f"by more than {P['spin tolerance']} are dropped."
            )
        text += (
            f" The files go to '{P['directory']}', with a record per molecule and "
            "a coverage table of the elements and bonded pairs."
        )
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self):
        """Export the training points."""
        next_node = super().run(printer)
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        printer.important(__(self.description_text(P), indent=self.indent))
        printer.important("")

        system_db = self.get_variable("_system_db")
        systems = self._systems(P, system_db)
        if not systems:
            raise ValueError("Export Training Set found no training systems.")

        directory = Path(P["directory"]).expanduser()
        if not directory.is_absolute():
            directory = Path(self.flowchart.root_directory) / directory
        directory.mkdir(parents=True, exist_ok=True)

        window = (
            None
            if self._is_none(P["energy window"])
            else self._as_float(P["energy window"], "kJ/mol")
        )
        spin = (
            None if self._is_none(P["spin tolerance"]) else float(P["spin tolerance"])
        )
        threshold = float(P["redraw threshold"]) / 100.0
        units = P["energy units"]

        coverage = {}
        summary = []
        refused = []
        redraw = []
        below = []  # molecules with points below their minimum
        n_kept = n_dropped = 0
        for name, first, second in self._molecules(systems):
            if first is None:
                first, second = second, None
            points, model, extra = self._points(P, first)
            redrawn = None
            if second is not None:
                r_points, r_model, _ = self._points(P, second)
                points, problem = export.merge_redraw(
                    points, r_points, same_model=(r_model == model)
                )
                if problem is not None:
                    refused.append((name, f"the redraw {second.name}: {problem}"))
                    continue
                for p in points:
                    if p.get("replaced"):
                        self._mark(p["configuration"], "replaced by the redraw", None)
                points = [p for p in points if not p.get("replaced")]
                redrawn = second.name
            kept, dropped, problem = export.filter_points(points, window, spin)
            if problem is not None:
                refused.append((name, problem))
                continue
            for p in kept:
                self._mark(p["configuration"], "kept", p["delta energy"])
            for (pname, kind, reason), conf in zip(
                dropped, [self._by_name(points, d[0]) for d in dropped]
            ):
                self._mark(conf, f"dropped: {reason}", None)
            fraction = export.sample_drop_fraction(points, dropped)
            if fraction > threshold:
                redraw.append((name, fraction))
            stem = self._file_stem(name)
            text = export.molecule_file(
                extra["symbols"],
                kept,
                name,
                units=units,
                charge=extra["charge"],
                multiplicity=extra["multiplicity"],
            )
            (directory / f"{stem}.extxyz").write_text(text)
            hessian_file = self._copy_hessian(extra["hessian"], directory, stem)
            record = export.provenance(
                name,
                model,
                kept[0].get("fingerprint"),
                extra["calculation"],
                extra["hessian"],
                extra["seed"],
                kept,
                dropped,
            )
            # the points refused while generating, by kind and reason
            record["dropped while generating"] = extra.get("generation drops") or {}
            if hessian_file is not None:
                record["hessian file"] = hessian_file
            record["sample drop fraction"] = round(fraction, 4)
            if redrawn is not None:
                record["samples from the redraw"] = redrawn
            lowest = min(p["delta energy"] for p in kept)
            record["lowest delta energy (kJ/mol)"] = round(lowest, 4)
            if lowest < -1.0:
                below.append((name, lowest))
            (directory / f"{stem}.json").write_text(export.dumps(record))
            bonds = sorted(generate.connectivity(extra["symbols"], kept[0]["xyz"]))
            export.add_coverage(
                coverage, extra["symbols"], bonds, [p["xyz"] for p in kept]
            )
            n_kept += len(kept)
            n_dropped += len(dropped)
            summary.append((name, len(kept), len(dropped), fraction))

        rows = export.coverage_table(
            coverage,
            (int(P["element minimum molecules"]), int(P["element minimum points"])),
            (int(P["pair minimum molecules"]), int(P["pair minimum points"])),
        )
        with open(directory / "coverage.csv", "w", newline="") as fd:
            writer = csv.writer(fd)
            writer.writerow(
                ["kind", "name", "molecules", "points", "r_min", "r_max", "passes"]
            )
            for row in rows:
                writer.writerow(
                    [
                        row[0],
                        row[1],
                        row[2],
                        row[3],
                        "" if row[4] is None else f"{row[4]:.4f}",
                        "" if row[5] is None else f"{row[5]:.4f}",
                        "yes" if row[6] else "no",
                    ]
                )

        data = {
            "number of molecules": len(summary),
            "number of points": n_kept,
            "number of dropped points": n_dropped,
            "refused molecules": [f"{n}: {r}" for n, r in refused],
            "molecules to redraw": [n for n, _ in redraw],
            # for a Loop over the systems to redraw ('name regexp'); matches nothing
            # when there are none
            "redraw regexp": (
                "^(?:" + "|".join(re.escape(n) for n, _ in redraw) + ")$"
                if redraw
                else "(?!)"
            ),
            "molecules below their minimum": [n for n, _ in below],
            "directory": str(directory),
        }
        self.analyze(
            data=data,
            summary=summary,
            refused=refused,
            redraw=redraw,
            below=below,
            rows=rows,
            threshold=threshold,
        )
        self.store_results(data=data)
        return next_node

    # ------------------------------------------------------------------ #
    def _systems(self, P, system_db):
        """The systems to export."""
        how = str(P["systems"]).strip()
        if how == "current":
            return [system_db.system]
        systems = list(system_db.systems)
        if how == "all training systems":
            return [s for s in systems if self._is_training_system(s)]
        return [s for s in systems if fnmatch.fnmatchcase(s.name, how)]

    @staticmethod
    def _molecules(systems):
        """(name, first pass, redraw) per molecule: a system '<x> redraw' is the
        second pass of samples for '<x> training'. Either may be None."""
        groups = {}
        order = []
        for system in systems:
            name = system.name
            if name.endswith(" redraw"):
                base, kind = name[: -len(" redraw")], "redraw"
            elif name.endswith(" training"):
                base, kind = name[: -len(" training")], "first"
            else:
                base, kind = name, "first"
            if base not in groups:
                groups[base] = {"first": None, "redraw": None, "name": None}
                order.append(base)
            groups[base][kind] = system
            if kind == "first":
                groups[base]["name"] = name
        out = []
        for base in order:
            g = groups[base]
            out.append((g["name"] or f"{base} training", g["first"], g["redraw"]))
        return out

    @staticmethod
    def _is_training_system(system):
        try:
            conf = system.configurations[0]
        except IndexError:
            return False
        return conf.properties.exists("point type#TrainingStructures#generation") and (
            _get(conf, "point type#TrainingStructures#generation") is not None
        )

    def _points(self, P, system):
        """The points of one training system, with their labels, and the
        molecule's facts (symbols, charge, multiplicity, records)."""
        configurations = list(system.configurations)
        model = P["model chemistry"]
        if model == "from the labels":
            models = set()
            for conf in configurations:
                for key in conf.properties.get().keys():
                    if key.startswith("energy#Energy#"):
                        models.add(key.split("#", 2)[2])
            if len(models) > 1:
                raise ValueError(
                    f"The points of {system.name} carry energies from "
                    f"{len(models)} model chemistries ({', '.join(sorted(models))}); "
                    "choose one in 'Labels from'."
                )
            model = models.pop() if models else None
        points = []
        for conf in configurations:
            points.append(
                {
                    "configuration": conf,
                    "name": conf.name,
                    "type": _get(conf, "point type#TrainingStructures#generation"),
                    "energy": _get(conf, f"energy#Energy#{model}"),
                    "gradients": _get(conf, f"gradients#Energy#{model}"),
                    "fingerprint": _get(conf, f"settings fingerprint#Energy#{model}"),
                    "S**2": _get(conf, f"S**2#Energy#{model}"),
                    "xyz": [
                        list(map(float, r))
                        for r in conf.atoms.get_coordinates(fractionals=False)
                    ],
                    "target energy": _get(
                        conf, "target energy#TrainingStructures#generation"
                    ),
                    "harmonic energy": _get(
                        conf, "harmonic energy#TrainingStructures#generation"
                    ),
                    "prescreen energy": _get(
                        conf, "prescreen energy#TrainingStructures#generation"
                    ),
                    "scan angle": _get(
                        conf, "scan angle#TrainingStructures#generation"
                    ),
                    "bond class": _get(
                        conf, "bond class#TrainingStructures#generation"
                    ),
                    "bond change": _get(
                        conf, "bond change#TrainingStructures#generation"
                    ),
                }
            )
        minimum = next(
            (p["configuration"] for p in points if p["type"] == "minimum"),
            configurations[0],
        )
        extra = {
            "symbols": list(minimum.atoms.symbols),
            "charge": int(minimum.charge or 0),
            "multiplicity": int(minimum.spin_multiplicity or 1),
            "calculation": _get(minimum, f"calculation#Energy#{model}"),
            "hessian": _get(minimum, "hessian#TrainingStructures#generation"),
            "seed": _get(minimum, "random seed#TrainingStructures#generation"),
            "generation drops": _get(
                minimum, "dropped points#TrainingStructures#generation"
            ),
        }
        return points, model, extra

    def _copy_hessian(self, hessian, directory, stem):
        """Copy the molecule's Hessian file (Training Structures' hessian.npz,
        whose path relative to the job is in its "hessian" record) next to its
        .extxyz as <stem>.hessian.npz; its name, or None if there is none."""
        import shutil

        name = (hessian or {}).get("file") if isinstance(hessian, dict) else None
        if not name:
            return None
        source = Path(name)
        if not source.is_absolute():
            source = Path(self.flowchart.root_directory) / source
        if not source.exists():
            printer.important(
                __(f"    Warning: the Hessian file {name} is missing.", indent="")
            )
            return None
        target = Path(directory) / f"{stem}.hessian.npz"
        shutil.copyfile(source, target)
        return target.name

    @staticmethod
    def _by_name(points, name):
        return next(p["configuration"] for p in points if p["name"] == name)

    def _mark(self, conf, status, delta):
        """Record on the configuration whether it was kept (never deleted)."""
        for name, value, _type, units in (
            ("training status#TrainingStructures#export", status, "str", None),
            ("delta energy#TrainingStructures#export", delta, "float", "kJ/mol"),
        ):
            if value is None:
                continue
            properties = conf.properties
            if not properties.exists(name):
                properties.add(name, _type, units=units, noerror=True)
            properties.put(name, value)

    @staticmethod
    def _file_stem(name):
        """A file name from a system name: runs of other characters -> '_'."""
        stem = re.sub(r"[^A-Za-z0-9._+@,-]+", "_", name).strip("_")
        return stem or "molecule"

    def analyze(
        self,
        data=None,
        summary=(),
        refused=(),
        redraw=(),
        below=(),
        rows=(),
        threshold=0.15,
        **kwargs,
    ):
        """Summarize the export."""
        if data is None:
            return
        printer.important(
            __(
                f"Exported {data['number of points']} points of "
                f"{data['number of molecules']} molecules to {data['directory']}; "
                f"{data['number of dropped points']} points were dropped.",
                indent=4 * " ",
            )
        )
        for name, problem in refused:
            printer.important(
                __(f"{name} was not exported: {problem}.", indent=4 * " ")
            )
        if redraw:
            names = ", ".join(f"{n} ({f:.0%})" for n, f in redraw)
            printer.important(
                __(
                    f"These molecules lost more than {threshold:.0%} of their samples "
                    "and need new samples with a lower maximum target energy (7 "
                    f"kcal/mol in the E2 recipe): {names}.",
                    indent=4 * " ",
                )
            )
        if below:
            names = ", ".join(f"{n} ({e:.1f} kJ/mol)" for n, e in below)
            printer.important(
                __(
                    "Some points lie below their molecule's minimum, which may not be "
                    f"the lowest conformer: {names}.",
                    indent=4 * " ",
                )
            )
        if rows:
            printer.important("")
            printer.important(__("Coverage of the exported points:", indent=4 * " "))
            lines = ["    kind     name     molecules   points   r range (Å)   passes"]
            for kind, name, m, n, r0, r1, ok in rows:
                r = "" if r0 is None else f"{r0:.3f}-{r1:.3f}"
                lines.append(
                    f"    {kind:8s} {name:8s} {m:9d} {n:8d}   {r:12s}  "
                    + ("yes" if ok else "no")
                )
            printer.important(__("\n".join(lines), indent=4 * " ", wrap=False))
        printer.important("")

    @staticmethod
    def _is_none(value):
        return isinstance(value, str) and value.strip().lower() in ("none", "")

    @staticmethod
    def _as_float(value, units):
        try:
            return value.m_as(units)
        except AttributeError:
            return float(value)


def _get(conf, name):
    """The value of a configuration's property, or None. The match is exact: names
    such as 'S**2#Energy#...' hold glob characters."""
    try:
        if not conf.properties.exists(name):
            return None
        found = conf.properties.get(name, match="exact")
    except Exception:
        return None
    entry = found.get(name) if isinstance(found, dict) else None
    if isinstance(entry, dict):
        return entry.get("value")
    return entry
