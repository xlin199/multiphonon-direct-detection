## Files

StructureFactor1D.py --
  1. Script calculating $S(\omega, q)$ as a function of $\omega$ for different $q$.
  2. Finished the incoherent approximation (see NSiteLatticeFixqRecursive.ipynb)

NSiteLatticeFixqRecursive.ipynb --
  1. Here focusing on plotting $S(\omega, q)$ as a function of $\omega$ for different $q$, using the recursive relation and incoherent approximation.

NSiteLatticeFixqRecursiveCoherent.ipynb --
  1. $S(\omega, q)$ as a function of $\omega$ for different $q$, including the off-diagonal factors
  2. Sanity check passed: 1 phonon resonance matches the dispersion relation

NSiteLatticeCoherentMultiprocess.py --
  1. Batch calculation of $S(\omega, q)$ for a list of $q$, same code as in NSiteLatticeFixqRecursiveCoherent.ipynb
   
  ---

/archive/NSiteLatticeFixq.ipynb --
  1. Here focusing on plotting $S(\omega, q)$ as a function of $\omega$ for different $q$
  2. --> ~~Why doesn't it match with the impulse approximation?~~ used wrong DW factor. Fixed in newer notebooks. 
  3. Moved to the /archive folder

NSiteLatticeFixq_2phonons.ipynb
  1. Visualizes the phonon energy distribution for the 2-phonon structure factor, incoherent corrections included
  2. 2d plot ( $\omega_1$ vs $\omega_2$ ) with colorbar $S^{(2)}(\omega=\omega_1+\omega_2, q)$
  3. Matches the dispersion relation

NSiteLattice1D.ipynb -- 
  1. Here focusing on plotting $\int S(\omega, q) d\omega$
  2. Directly calculating the structure factor without any smearing of the phonon 
  distribution, and without using the recursion relation.
  3. --> Plots comparing the incoherent approx with the full coherent result - still useful! 
  4. --> Can do one version with energy cut and one without energy cut (N-dependent DW)
  
NSiteLattice2D.ipynb -- 
  1. Here focusing on plotting $S(\omega, q)$
  2. --> Still useful to have 2D plots of $S(\omega, q)$! [But first focus on the 1D plots]
  3. Currently old independent code from Momei -- would either need to perform checks against
  4. the code in NSiteLattice1D.ipynb, or can rewrite it given the latest results in calculating structure

---

<!-- TODO:
  1. Update notebooks for the N-site lattice 
       * $\int S(\omega, q) d\omega$ vs. q comparing incoherent and coherent
       * Plotting $S(\omega, q)$ as a function of $\omega$ for different $q$ 
          -- compare incoherent and coherent, compare impulse approx
       * 2D plot $S(\omega, q)$
       * Rate calculation (compare with single HO)
  . New notebook for reproducing single HO structure factor and rate
  1. Notebook for structure factor and rate for two-site lattice
       * $\int S(\omega, q) d\omega$ vs. q comparing incoherent and coherent
       * Plotting $S(\omega, q)$ as a function of $\omega$ for different $q$ 
          -- compare incoherent and coherent, compare impulse approx
       * 2D plot $S(\omega, q)$
       * Rate calculation (compare with single HO) -->

#### TODO
1. Rate calculation 
   - Sanity check: reproduce 2205.02250 Fig 15(a), silicon with 1 meV energy threshold (intermediate step: Fig 13 (c) $dR/d \omega$)


