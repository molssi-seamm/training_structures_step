.. _user-guide:

**********
User Guide
**********

..
    <remove the dots above and this line and unindent the toctree to expose it>
    Contents:

    .. toctree::
       :glob:
       :maxdepth: 2
       :titlesonly:

       *

What the step makes
===================

The step works on the current configuration, which should be a minimum of the model
chemistry chosen by a *Model Chemistry* step earlier in the flowchart. From the Hessian
of that model chemistry it makes, in a new system named '<molecule> training':

- **the minimum** itself, as the first configuration;
- **normal-mode samples** ('sample 1', ...): the classical shape, which gives every mode
  the same expected energy, scaled to a harmonic energy drawn uniformly between the
  minimum and maximum target energies (0.5 and 9 kcal/mol by default). Modes below the
  lowest frequency (100 cm⁻¹) are left to the scans. By default there are half as many
  samples as sampled modes, between 10 and 40. These are the samples of the Normal Mode
  Sampling step's energy-targeted distribution;
- **relaxed torsion scans** ('scan 1 C2-C3 +45', ...) of up to two heavy rotors: acyclic,
  non-conjugated single bonds with a heavy atom beyond each end, those that move the most
  heavy atoms first. Methyl, OH and NH₂ rotations are left to the samples, and conjugated
  bonds (amides, esters, aryl ketones) are left out because their barriers are high. The
  dihedral goes around the circle in 15° steps; each point is relaxed with the dihedral
  held fixed, at the generator level (GFN2-xTB by default);
- **bond stretches** ('stretch C-H +', 'stretch C-H -'): the median bond of each bonded
  element pair is lengthened and shortened by a, where ½ a² vᵀHv equals the stretch
  energy (8 kcal/mol), at most 0.10 Å;
- for a **diatomic**, a stretch grid of 100 points out to a harmonic energy of
  9.5 kcal/mol, in place of the samples and stretches.

Checks
======

Every point is checked before it is kept:

- **The bonds**: atoms closer than 1.3 times the sum of their covalent radii are bonded,
  and a point must have exactly the minimum's bonds (a stretch point's own bond
  excepted). A sample that fails is drawn again; other points are dropped.
- **The pre-screen**: the energy at the generator level must be within 15 kcal/mol of
  the minimum's.
- **The limits on samples**: no atom moves more than 0.30 Å, no bond length changes by
  more than 0.20 Å and no bond angle by more than 15°. A sample is drawn up to five
  times, keeping its target energy, then dropped.

The output lists the points dropped and why.

The Hessian
===========

For open-shell molecules, and a random 5 % of closed-shell ones, the Hessian is checked
along three random displacements δ (0.005 Å RMS) against central differences of the
gradient: e = ‖Hδ − [g(x₀+δ) − g(x₀−δ)]/2‖ / ‖Hδ‖. If any e exceeds 0.02, or any
frequency exceeds 5000 cm⁻¹, the Hessian is refused and replaced by finite differences
of the gradient. How the Hessian was made and checked is stored with the minimum.

Properties of each point
========================

Each configuration carries, as properties ``<name>#TrainingStructures#generation``:
the point type; the target and harmonic energies; the pre-screen energy; for scans the
dihedral and the change of angle; for stretches the element pair, the bond and the
change of length. The random seed is reported and saved with the results.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
