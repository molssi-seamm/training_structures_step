=======
History
=======

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
