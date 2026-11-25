# -*- coding: utf-8 -*-
"""
Control parameters for the Training Structures step in a SEAMM flowchart
"""

import logging
import seamm
import pprint  # noqa: F401

logger = logging.getLogger(__name__)


def _p(default, kind, units, description, help_text, enumeration=(), fmt=""):
    return {
        "default": default,
        "kind": kind,
        "default_units": units,
        "enumeration": tuple(enumeration),
        "format_string": fmt,
        "description": description,
        "help_text": help_text,
    }


class TrainingStructuresParameters(seamm.Parameters):
    """
    The control parameters for Training Structures.

    The step makes distorted structures of one molecule (the current
    configuration, which should be a minimum of the flowchart's Model Chemistry)
    for training a machine-learned potential: normal-mode samples at targeted
    energies, relaxed torsion scans of the heavy rotors, and bond stretches. The
    Hessian comes from the Model Chemistry; the scans and the pre-screen use a
    second, cheap level.
    """

    parameters = {
        # ------------------------------------------------------------ samples
        "number of samples": _p(
            "in proportion to the modes",
            "integer",
            "",
            "Number of samples:",
            "How many normal-mode samples to draw: a number, or 'in proportion to "
            "the modes' for the number of sampled modes times 'samples per mode', "
            "kept between the minimum and maximum below. 0 makes no samples.",
            ("in proportion to the modes",),
        ),
        "samples per mode": _p(
            0.5,
            "float",
            "",
            "Samples per mode:",
            "Samples per sampled mode.",
            fmt=".2f",
        ),
        "minimum number of samples": _p(
            10, "integer", "", "At least:", "The fewest samples."
        ),
        "maximum number of samples": _p(
            40, "integer", "", "At most:", "The most samples."
        ),
        "minimum target energy": _p(
            0.5,
            "float",
            "kcal/mol",
            "Target energies from:",
            "Each sample is the equipartition shape scaled to a harmonic energy "
            "drawn uniformly between this and the maximum.",
            fmt=".2f",
        ),
        "maximum target energy": _p(
            9.0, "float", "kcal/mol", "to:", "The highest target energy.", fmt=".2f"
        ),
        "amplitude cap": _p(
            3.0,
            "float",
            "",
            "Amplitude cap (in sigma):",
            "Clamp each normal coordinate of the shape to this many standard "
            "deviations before scaling. 'none' leaves the tails uncapped.",
            ("none",),
            ".1f",
        ),
        "lowest frequency": _p(
            100.0,
            "float",
            "1/cm",
            "Skip modes below:",
            "Do not sample along modes softer than this: soft torsions are covered "
            "by the scans. 'none' samples every mode.",
            ("none",),
            ".1f",
        ),
        "maximum atom displacement": _p(
            0.30,
            "float",
            "Å",
            "Largest atom displacement:",
            "Redraw a sample in which an atom moves farther than this. 'none' sets "
            "no limit.",
            ("none",),
            ".3f",
        ),
        "maximum bond change": _p(
            0.20,
            "float",
            "Å",
            "Largest bond-length change:",
            "Redraw a sample in which a bond length changes by more than this.",
            ("none",),
            ".3f",
        ),
        "maximum angle change": _p(
            15.0,
            "float",
            "degree",
            "Largest bond-angle change:",
            "Redraw a sample in which a bond angle changes by more than this.",
            ("none",),
            ".1f",
        ),
        "maximum draws": _p(
            5,
            "integer",
            "",
            "Draws per sample:",
            "How many times a sample is drawn (keeping its target energy) before "
            "it is dropped, when it fails the limits, the connectivity check or the "
            "pre-screen.",
        ),
        # ------------------------------------------------------------ scans
        "number of rotors": _p(
            2,
            "integer",
            "",
            "Rotors to scan:",
            "Scan up to this many heavy rotors (acyclic, non-conjugated single "
            "bonds with a heavy atom beyond each end), those that move the most "
            "heavy atoms first. 0 makes no scans.",
        ),
        "scan step": _p(
            15.0,
            "float",
            "degree",
            "Scan step:",
            "The step of each torsion scan; the minimum's own angle is not "
            "repeated.",
            fmt=".1f",
        ),
        "generator level": _p(
            "xTB:SQM@GFN2-xTB",
            "string",
            "",
            "Level for scans and pre-screen:",
            "The cheap model chemistry that relaxes the scan points at fixed "
            "torsion and pre-screens every point. It only makes geometries: no "
            "energy of it is a label.",
            ("xTB:SQM@GFN2-xTB", "MOPAC:SQM@PM7"),
        ),
        "relaxation convergence": _p(
            0.01,
            "float",
            "eV/Å",
            "Relax scan points to:",
            "The largest force on an atom (the torsion held fixed) at which a scan "
            "point's relaxation stops.",
            fmt=".3f",
        ),
        "relaxation max steps": _p(
            200, "integer", "", "At most:", "Steps of each scan point's relaxation."
        ),
        # ------------------------------------------------------------ bonds
        "bond stretches": _p(
            "yes",
            "enum",
            "",
            "Bond stretches:",
            "Stretch and compress the median bond of each bonded element pair by "
            "the amplitude below (two points per pair).",
            ("yes", "no"),
        ),
        "bond pairs": _p(
            "all",
            "string",
            "",
            "Element pairs:",
            "'all', or the pairs to stretch, e.g. 'C-H, C-C, C-O'.",
            ("all",),
        ),
        "stretch energy": _p(
            8.0,
            "float",
            "kcal/mol",
            "Stretch to an energy of:",
            "The amplitude a solves 1/2 a² vᵀHv = this, with the full Hessian.",
            fmt=".2f",
        ),
        "maximum stretch": _p(
            0.10,
            "float",
            "Å",
            "but by at most:",
            "The largest bond-length change of a stretch point.",
            fmt=".3f",
        ),
        "diatomic grid points": _p(
            100,
            "integer",
            "",
            "Diatomic stretch grid:",
            "A diatomic gets this many points evenly spaced along its stretch, in "
            "place of samples and stretches.",
        ),
        "diatomic grid energy": _p(
            9.5,
            "float",
            "kcal/mol",
            "out to a harmonic energy of:",
            "The ends of a diatomic's grid.",
            fmt=".2f",
        ),
        # ------------------------------------------------------------ filters
        "connectivity check": _p(
            "yes",
            "enum",
            "",
            "Keep the bonds:",
            "Refuse a point whose bonds (atoms closer than the factor below times "
            "the sum of their covalent radii) differ from the minimum's.",
            ("yes", "no"),
        ),
        "connectivity factor": _p(
            1.3, "float", "", "Bond factor:", "Covalent-radius factor.", fmt=".2f"
        ),
        "prescreen energy": _p(
            15.0,
            "float",
            "kcal/mol",
            "Pre-screen energy:",
            "Refuse a point whose energy at the generator level is more than this "
            "above the minimum's. 'none' skips the pre-screen.",
            ("none",),
            ".1f",
        ),
        # ------------------------------------------------------------ Hessian
        "check the hessian": _p(
            "open shells",
            "enum",
            "",
            "Check the Hessian:",
            "Test the Hessian against central differences of the gradient along "
            "random displacements; a refused one is replaced by finite "
            "differences of the gradient. 'open shells' checks every open-shell "
            "molecule and the fraction below of closed-shell ones.",
            ("open shells", "always", "never"),
        ),
        "closed-shell check fraction": _p(
            0.05,
            "float",
            "",
            "and closed shells at random:",
            "The fraction of closed-shell molecules checked, as a watch.",
            fmt=".2f",
        ),
        "hessian tolerance": _p(
            0.02,
            "float",
            "",
            "Relative error allowed:",
            "e = |Hδ - Δg/2| / |Hδ| above this, or any |ν| above the maximum "
            "frequency, refuses the Hessian.",
            fmt=".3f",
        ),
        "check displacement": _p(
            0.005,
            "float",
            "Å",
            "Check displacement (RMS):",
            "The size of each random displacement.",
            fmt=".4f",
        ),
        "maximum frequency": _p(
            5000.0,
            "float",
            "1/cm",
            "Maximum frequency:",
            "A frequency above this means the Hessian is broken.",
            fmt=".0f",
        ),
        # ------------------------------------------------------------ output
        "random seed": _p(
            "random",
            "string",
            "",
            "Random seed:",
            "An integer for reproducible points, or 'random' (the seed used is "
            "recorded).",
            ("random",),
        ),
        "system name": _p(
            "from structure",
            "string",
            "",
            "Name of the new system:",
            "The system holding the points: 'from structure' adds ' training' to "
            "the molecule's name.",
            ("from structure",),
        ),
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("TrainingStructuresParameters.__init__")

        super().__init__(
            defaults={**TrainingStructuresParameters.parameters, **defaults}, data=data
        )
