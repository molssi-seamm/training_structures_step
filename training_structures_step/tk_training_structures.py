# -*- coding: utf-8 -*-

"""The graphical part of a Training Structures step"""

import pprint  # noqa: F401
import tkinter as tk  # noqa: F401
import tkinter.ttk as ttk

import training_structures_step  # noqa: F401
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw


class TkTrainingStructures(seamm.TkNode):
    """
    The graphical part of a Training Structures step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    namespace : str
        The namespace of the current step.
    tk_subflowchart : TkFlowchart
        A graphical Flowchart representing a subflowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in Training Structures_parameters.py

    See Also
    --------
    TrainingStructures, TkTrainingStructures,
    TrainingStructuresParameters,
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        namespace: str
            The stevedore namespace for finding sub-nodes.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )

    #: the dialog's frames and the parameters in each, in order
    GROUPS = {
        "Samples": (
            "number of samples",
            "samples per mode",
            "minimum number of samples",
            "maximum number of samples",
            "minimum target energy",
            "maximum target energy",
            "amplitude cap",
            "lowest frequency",
            "maximum atom displacement",
            "maximum bond change",
            "maximum angle change",
            "maximum draws",
        ),
        "Torsion scans": (
            "number of rotors",
            "scan step",
            "relaxation convergence",
            "relaxation max steps",
        ),
        "Bond stretches": (
            "bond stretches",
            "bond pairs",
            "stretch energy",
            "maximum stretch",
            "diatomic grid points",
            "diatomic grid energy",
        ),
        "Checks": (
            "generator level",
            "connectivity check",
            "connectivity factor",
            "prescreen energy",
        ),
        "Hessian": (
            "check the hessian",
            "closed-shell check fraction",
            "hessian tolerance",
            "check displacement",
            "maximum frequency",
        ),
        "Output": ("random seed", "system name"),
    }
    #: the parameters whose value changes the layout
    DRIVERS = (
        "number of samples",
        "number of rotors",
        "bond stretches",
        "connectivity check",
        "prescreen energy",
        "check the hessian",
    )

    def create_dialog(self):
        """Create the dialog: the parameters in labelled frames."""
        frame = super().create_dialog(title="Training Structures")
        P = self.node.parameters
        for group, keys in self.GROUPS.items():
            self[group] = ttk.LabelFrame(
                frame, borderwidth=4, relief="sunken", text=group, labelanchor="n"
            )
            for key in keys:
                self[key] = P[key].widget(self[group])
        self["model chemistry note"] = ttk.Label(
            frame,
            text="Add a Model Chemistry step before this one: it gives the Hessian.",
            foreground="red",
        )
        for key in self.DRIVERS:
            w = self[key]
            entry = w.combobox if hasattr(w, "combobox") else w.entry
            for event in ("<<ComboboxSelected>>", "<Return>", "<FocusOut>"):
                entry.bind(event, self.reset_dialog)
        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Lay out the frames, showing each parameter only when it applies."""
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        def value(key):
            return str(self[key].get()).strip()

        shown = {key for keys in self.GROUPS.values() for key in keys}
        if value("number of samples") != "in proportion to the modes":
            shown -= {
                "samples per mode",
                "minimum number of samples",
                "maximum number of samples",
            }
        if value("number of samples") == "0":
            shown -= set(self.GROUPS["Samples"][1:])
        if value("number of rotors") == "0":
            shown -= set(self.GROUPS["Torsion scans"][1:])
        if value("bond stretches") == "no":
            shown -= {"bond pairs", "stretch energy", "maximum stretch"}
        if value("connectivity check") == "no":
            shown.discard("connectivity factor")
        if value("number of rotors") == "0" and value("prescreen energy") == "none":
            shown.discard("generator level")
        if value("check the hessian") == "never":
            shown -= set(self.GROUPS["Hessian"][1:])
        elif value("check the hessian") == "always":
            shown.discard("closed-shell check fraction")

        row = 0
        if not self._upstream_has_model_chemistry():
            self["model chemistry note"].grid(row=row, column=0, sticky=tk.W)
            row += 1
        for group, keys in self.GROUPS.items():
            inner = self[group]
            for slave in inner.grid_slaves():
                slave.grid_forget()
            widgets = []
            for key in keys:
                if key in shown:
                    self[key].grid(row=len(widgets), column=0, sticky=tk.EW)
                    widgets.append(self[key])
            if widgets:
                sw.align_labels(widgets, sticky=tk.E)
                inner.grid(row=row, column=0, sticky=tk.EW, pady=4)
                row += 1

        have_results = (
            "results" in self.node.metadata and len(self.node.metadata["results"]) > 0
        )
        if have_results and "results" in self.node.parameters:
            self.setup_results()

    def _upstream_has_model_chemistry(self):
        """True if a Model Chemistry step precedes this one (False on any
        error, so the reminder shows)."""
        try:
            return any(
                type(node).__name__ == "ModelChemistry"
                and type(node).__module__.startswith("model_chemistry_step")
                for node in self.previous_nodes()
            )
        except Exception:
            return False

    def right_click(self, event):
        """
        Handles the right click event on the node.

        Parameters
        ----------
        event : Tk Event

        Returns
        -------
        None

        See Also
        --------
        TkTrainingStructures.edit
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
