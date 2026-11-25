=======
History
=======

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
