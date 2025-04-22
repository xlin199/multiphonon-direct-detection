NSiteLattice1D.ipynb -- 
  Here focusing on plotting $\int S(\omega, q) d\omega$
  Directly calculating the structure factor without any smearing of the phonon 
  distribution, and without using the recursion relation.
  --> Plots comparing the incoherent approx with the full coherent result - still useful! 
  --> Can do one version with energy cut and one without energy cut (N-dependent DW)
  
NSiteLattice2D.ipynb -- 
  Here focusing on plotting $S(\omega, q)$
  --> Still useful to have 2D plots of $S(\omega, q)$! [But first focus on the 1D plots]
  Currently old independent code from Momei -- would either need to perform checks against
  the code in NSiteLattice1D.ipynb, or can rewrite it given the latest results in calculating structure

NSiteLatticeFixqRecursive.ipynb --
  Here focusing on plotting $S(\omega, q)$ as a function of $\omega$ for different $q$
  --> Why doesn't it match with the impulse approximation?
  [Can move NSiteLatticeFixq.ipynb to a folder of old notebooks]

TO DO:
  1. Move "class StructureFactor1D" to separate file, have all the notebooks just calling this and making the various plots
  2. Inside "class StructureFactor1D", we can have the recursive structure factor calculation 
      and also keep the direct calculation in NSiteLattice1D.ipynb
  3. Update notebooks for the N-site lattice 
       * $\int S(\omega, q) d\omega$ vs. q comparing incoherent and coherent
       * Plotting $S(\omega, q)$ as a function of $\omega$ for different $q$ 
          -- compare incoherent and coherent, compare impulse approx
       * 2D plot $S(\omega, q)$
       * Rate calculation (compare with single HO)
  4. New notebook for reproducing single HO structure factor and rate
  5. Notebook for structure factor and rate for two-site lattice
       * $\int S(\omega, q) d\omega$ vs. q comparing incoherent and coherent
       * Plotting $S(\omega, q)$ as a function of $\omega$ for different $q$ 
          -- compare incoherent and coherent, compare impulse approx
       * 2D plot $S(\omega, q)$
       * Rate calculation (compare with single HO)
