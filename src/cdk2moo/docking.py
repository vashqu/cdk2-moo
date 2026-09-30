"""Stage 6: docking with AutoDock Vina. ISOLATED ON PURPOSE.

    Nothing outside this module and the docking scripts may import from here.
    Stage 6 is timeboxed and droppable: if docking is not working by the Day 3
    midpoint it gets deleted entirely, and stages 1-5 and 8 must continue to
    run untouched. Keep the dependency arrow pointing one way.

Will live here:

- receptor preparation via Meeko's ``mk_prepare_receptor.py --read_pdb``
  (never ADFR Suite or MGLTools: no arm64 builds)
- ligand preparation to PDBQT, including the N-oxide in dinaciclib, which is
  a known atom-typing trap
- box definition centred on the co-crystallised ligand
- Vina invocation, with replicate runs under different seeds so engine noise
  is measurable
- redocking of the native ligand and RMSD to the crystal pose, which is run
  FIRST; above 2 A, every downstream docking number is reported as untrusted
"""
