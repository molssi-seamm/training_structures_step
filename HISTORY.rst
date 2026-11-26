=======
History
=======

2026.11.26.3 -- Bugfix: no rotors onto linear groups
    * A single bond onto a linear centre (a nitrile's or alkyne's sp carbon, an
      isocyanide's nitrogen, an allene's central carbon, or atoms that are collinear
      in the structure) is no longer a rotor. Its torsion is undefined, and scanning
      it stopped the molecule with "FixInternals.adjust_positions did not converge"
      (ethyl cyanide, 2-cyclopentenecarbonitrile, 1-norbornylisocyanide).
    * A scan point whose torsion cannot be held is dropped and reported, rather than
      stopping the molecule.
    * The points refused while generating are stored on the minimum ("dropped
      points", by kind and reason, with counts) and the export's per-molecule record
      carries them as "dropped while generating".

2026.11.26.2 -- A replacement Hessian is checked too
    * The Hessian check's gradients each start from the reference structure's
      orbitals, with a tight SCF (seamm-exec 2026.11.26).
    * A refused Hessian's finite-difference replacement, made the same way with every
      point checked for a state change, is itself checked; if it fails too the
      molecule is refused. Before, the replacement was used unchecked (methoxy in the
      E2b pilot, SEAMM job 6090).
    * How the finite differences were made, both checks and the zero modes are in the
      molecule's record.

2026.11.26.1 -- Record what the engine ran for the Hessian
    * Training Structures stores, with how the Hessian was made, the warm engine's own
      record of the Hessian run (else its last energy run): the settings fingerprint, the
      program version and the input's '!' line(s). The engine runs in node-local scratch,
      so its input and output are not kept; this record, made from them, is what shows
      the Hessian's level and settings afterwards. Export Training Set carries it in each
      molecule's record.

2026.11.26 -- Choose how the Hessian is made; the warm engine by default
    * A new parameter, "Energies and Hessian via:", as the Thermochemistry step's: the
      warm MDI engine (analytic where the program has it, and needed for the Hessian
      check), separate tasks (central finite differences of 6N gradients where the job
      runs its tasks), or automatic.
    * The default is the warm engine, as the Thermochemistry step's CBH approach. On a
      queue that runs tasks, the old automatic choice made every Hessian 6N gradients as
      tasks, several times the cost of the molecule's labels, and the Hessian check
      could not run.
    * The dialog hides the check's settings for separate tasks, which are already finite
      differences.

2026.11.25.2 -- A redraw pass for molecules that lost too many samples
    * Training Structures can name its system 'for a redraw' ('<molecule> redraw'),
      for a second pass of samples from the first pass's minimum, e.g. with a lower
      maximum target energy, no scans and no stretches.
    * Export Training Set merges a molecule's redraw with its first pass: the redraw's
      samples replace the first pass's (marked 'replaced by the redraw'), the scans and
      stretches are kept, and both passes must share the minimum's energy and settings.
      The record says where the samples came from.
    * Export Training Set gives the molecules to redraw as a regular expression
      ('redraw regexp'), for a Loop over those systems in the same flowchart.

2026.11.25.1 -- Export Training Set: labelled points to extended-XYZ files
    * A new step, Export Training Set, for the points of Training Structures after an
      Energy step has labelled them: each point's energy above its molecule's minimum, a
      check that a molecule's points share one settings fingerprint (no mixed levels),
      and filters for an energy window (10 kcal/mol) and, for open shells, a change of
      <S**2> (0.05). Dropped points are marked on their configurations, not deleted.
    * One extended-XYZ file per molecule (REF_energy, REF_forces, delta_energy and the
      point's details, full-precision positions; keys that ASE reads as data, not as a
      calculator's results), with a JSON record of how it was made: model chemistry,
      settings fingerprint, the calculation record, how the Hessian was made and
      checked, the random seed, and the points kept and dropped.
    * A coverage table of elements and bonded pairs (molecules, points, the bond lengths
      sampled) against minimum counts; molecules that lost too many samples are
      reported for a redraw, and points below a molecule's minimum are reported.
    * Training Structures now stores its random seed on the minimum.

2026.11.25 -- Initial release
    * Distorted structures of a molecule for training a machine-learned potential, from
      the Hessian of the Model Chemistry: energy-targeted normal-mode samples, relaxed
      torsion scans of up to a set number of heavy rotors at a cheap generator level,
      and two stretches of the median bond of each bonded element pair, or a stretch
      grid for a diatomic.
    * Every point must keep the molecule's bonds and lie within a set energy of the
      minimum at the generator level; refused points are counted with the reasons.
    * The Hessian is checked against finite differences of the gradient for open-shell
      molecules (and a fraction of closed-shell ones) and replaced by a
      finite-difference Hessian if it fails.
    * The points are configurations of one new system, named by kind, with their
      generation details as properties.
