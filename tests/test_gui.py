# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialog: create it and re-lay it out for every choice
that drives the layout. Skipped when no display is available."""

import pytest


@pytest.fixture()
def tk_node():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw
    import seamm

    Pmw.initialise(root)
    flowchart = seamm.Flowchart(namespace="org.molssi.seamm", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.tk"
    )
    node = flowchart.create_node("Training Structures")
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get("Training Structures")
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    yield tk_node
    root.destroy()


def shown(tk_node, key):
    return tk_node[key].grid_info() != {}


def test_dialog_layouts(tk_node):
    tk_node.create_dialog()
    tk_node.reset_dialog()
    # defaults: in proportion to the modes, 2 rotors, stretches, checks on
    for key in (
        "samples per mode",
        "scan step",
        "stretch energy",
        "generator level",
        "connectivity factor",
        "closed-shell check fraction",
    ):
        assert shown(tk_node, key), key
    tk_node["number of samples"].set("20")
    tk_node.reset_dialog()
    assert not shown(tk_node, "samples per mode")
    assert shown(tk_node, "minimum target energy")
    tk_node["number of samples"].set("0")
    tk_node.reset_dialog()
    assert not shown(tk_node, "minimum target energy")

    tk_node["number of rotors"].set("0")
    tk_node.reset_dialog()
    assert not shown(tk_node, "scan step")
    assert shown(tk_node, "generator level")  # still needed for the pre-screen
    tk_node["prescreen energy"].set("none")
    tk_node.reset_dialog()
    assert not shown(tk_node, "generator level")

    tk_node["bond stretches"].set("no")
    tk_node.reset_dialog()
    assert not shown(tk_node, "stretch energy")
    assert shown(tk_node, "diatomic grid points")

    tk_node["connectivity check"].set("no")
    tk_node.reset_dialog()
    assert not shown(tk_node, "connectivity factor")

    tk_node["check the hessian"].set("never")
    tk_node.reset_dialog()
    assert not shown(tk_node, "hessian tolerance")
    tk_node["check the hessian"].set("always")
    tk_node.reset_dialog()
    assert shown(tk_node, "hessian tolerance")
    assert not shown(tk_node, "closed-shell check fraction")
