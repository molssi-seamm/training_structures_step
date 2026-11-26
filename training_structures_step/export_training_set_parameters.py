# -*- coding: utf-8 -*-
"""
Control parameters for the Export Training Set step in a SEAMM flowchart
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


class ExportTrainingSetParameters(seamm.Parameters):
    """
    The control parameters for Export Training Set.

    The step reads the systems made by Training Structures after an Energy step
    has labelled them, works out each point's energy above its molecule's
    minimum, drops points outside the energy window or with a changed <S**2>, and
    writes one extended-XYZ file per molecule, with a record of how it was made,
    and a coverage table of the elements and bonded pairs.
    """

    parameters = {
        "systems": _p(
            "all training systems",
            "string",
            "",
            "Systems:",
            "Which systems to export: 'all training systems' (every system whose "
            "configurations were made by Training Structures), 'current', or a "
            "pattern for the system names (e.g. '* training').",
            ("all training systems", "current"),
        ),
        "model chemistry": _p(
            "from the labels",
            "string",
            "",
            "Labels from:",
            "The model chemistry whose energies and gradients are the labels: "
            "'from the labels' when the points carry only one, else its level "
            "(e.g. 'ORCA:DFT@r2SCAN-3c').",
            ("from the labels",),
        ),
        "energy window": _p(
            10.0,
            "float",
            "kcal/mol",
            "Energy window:",
            "Drop a point more than this above its molecule's minimum. 'none' "
            "keeps every labelled point.",
            ("none",),
            ".2f",
        ),
        "spin tolerance": _p(
            0.05,
            "float",
            "",
            "Allowed change of <S**2>:",
            "For an open-shell molecule, drop a point whose <S**2> differs from the "
            "minimum's by more than this (a different electronic state or more spin "
            "contamination). 'none' skips the check.",
            ("none",),
            ".3f",
        ),
        "redraw threshold": _p(
            15.0,
            "float",
            "",
            "Report a redraw above (%):",
            "A molecule that lost more than this share of its samples to the "
            "window or the spin check is reported for new samples with a lower "
            "maximum target energy.",
            fmt=".1f",
        ),
        "energy units": _p(
            "eV",
            "enum",
            "",
            "Units in the file:",
            "The units of REF_energy and delta_energy (forces in the same per Å).",
            ("eV", "kJ/mol", "E_h"),
        ),
        "directory": _p(
            "training_set",
            "string",
            "",
            "Directory:",
            "Where the files go, relative to the job's directory unless absolute: "
            "<molecule>.extxyz and <molecule>.json per molecule, and coverage.csv.",
        ),
        "element minimum molecules": _p(
            30, "integer", "", "Elements need at least:", "molecules ..."
        ),
        "element minimum points": _p(500, "integer", "", "and points:", "points."),
        "pair minimum molecules": _p(
            10, "integer", "", "Bonded pairs need at least:", "molecules ..."
        ),
        "pair minimum points": _p(100, "integer", "", "and points:", "points."),
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
        logger.debug("ExportTrainingSetParameters.__init__")
        super().__init__(
            defaults={**ExportTrainingSetParameters.parameters, **defaults}, data=data
        )
