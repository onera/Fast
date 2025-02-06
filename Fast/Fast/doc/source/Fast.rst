.. Fast documentation master file

:tocdepth: 2

Fast: Flexible Aerodynamic Solver Technolgy
============================================

Preamble
########

Fast is a common/front module for all Fast series solvers.

Fast is only available for use with the pyTree interface. You must 
import the module::

    import Fast.PyTree as Fast


.. py:module:: Fast

List of functions
##################

**-- Actions**

.. autosummary::
   :nosignatures:
   
   Fast.PyTree.setNum2Base
   Fast.PyTree.setNum2Zones
   Fast.PyTree.loadTree
   Fast.PyTree.load
   Fast.PyTree.saveTree
   Fast.PyTree.save
   
.. Fast.metric ..



Contents
#########

Actions
--------------------------

.. py:function:: Fast.PyTree.setNum2Base(a, numb)

    Set the numb dictionary to bases. Exists also
    as in place version (_setNum2Base) that modifies a
    and returns None.

    :param a: input data
    :type a: Base, pyTree
    :param numb: numerics for base
    :type numb: dictionary

    the **keys** of **numb** dictionary are:

    - **'temporal_scheme'**: numerical scheme for time integration. Possible values are

        + 'explicit' (RK3 scheme, see p49 http://publications.onera.fr/exl-php/util/documents/accede_document.php)
        + 'implicit' (BDF2 (second order Gear scheme for global time stepping) or BDF1 (first order Euler scheme for local time stepping))
        + 'implicit_local' (see p107 http://publications.onera.fr/exl-php/docs/ILS_DOC/227155/DOC356618_s1.pdf)
        + default value is 'implicit'

    - **'ss_iteration'**: Newton iterations for implicit schemes.

        + default value is 30

    - **'modulo_verif'**: computation period for cfl (RK3 or BDF2), Newton convergence (all), and predictor estimation (implicit_local only)

        + default value is 200

    *Example of use:*

    * `Set numerics to base (pyTree) <Examples/Fast/setNum2BasePT.py>`_:

    .. literalinclude:: ../../test/setNum2BasePT.py

------------------------------------------------------------------

.. py:function:: Fast.PyTree.setNum2Zones(a, numz)

    Set the numz dictionary to zones. Exists also
    as in place version (_setNum2Zones) that modifies a
    and returns None.
    
    :param a: input data
    :type a: Zone, Base, pyTree
    :param numz: input data
    :type numz: dictionary

    the **keys** of **numz** dictionary are:

    - **'scheme'**: numerical scheme for convective flux reconstruction. Possible values are
 
        + 'ausmpred' (AUSM+(P) reduced scheme, see p49 https://tel.archives-ouvertes.fr/pastel-00834850/document)
        + 'senseur'  (AUSM+(P)-based hybrid centered/decentered sensor scheme for DNS/LES, see p50 https://tel.archives-ouvertes.fr/pastel-00834850/document)
        + 'roe'      (classical Roe scheme)
        + default value is 'ausmpred'

    - **'slope'**: Slope reconstruction method. Possible values are

        + 'o1'      (first order MUSCL reconstrution, only valid for roe scheme)
        + 'minmod'  (second order MUSCL reconstrution with minmod limiter, only valid for roe scheme)
        + 'o3'      (third order MUSCL reconstrution, see  p50 https://tel.archives-ouvertes.fr/pastel-00834850/document)
        + 'o3sc'    (third order MUSCL reconstrution with slope limiter, see :download:`Lugrin scheme <documents/SCHEMA_lugrin.pdf>`)
        + 'o5'      (fifth order MUSCL reconstrution)
        + 'o5sc'    (fifth order MUSCL reconstrution with slope limiter)
        + default value is 'o3'

    - **'senseurType'**:  Sensor mode for the 'senseur' scheme. Possible values are

        + 0  (correction for speed only)
        + 1  (correction for speed, density and pressure)
        
    - **'motion'**: Motion mode for moving walls. Possible values are

        + 'none'        (no motion)
        + 'rigid'       (rigid motion defined in Fast, see p47 https://tel.archives-ouvertes.fr/tel-01011273/document)
        + 'rigid_ext'   (rigid motion defined externally by RigidMotion)
        + 'deformation' (ALE with deformation)
        + default value is 'none'

    - **'time_step'**: time step value

        + default value is 1e-4

    - **'time_step_nature'**: time step nature. Possible values are 

        + 'global'
        + 'local'
        + default value is 'global'

    - **'cfl'**: Courant-Friedrichs-Lewy condition (only for time_step_nature'='local')

        + default value is 1

    - **'epsi_newton'**: Newton stopping criteria based on Loo norm

        + default value is 0.1
    
    - **'inj1_newton_tol'**: Newton tolerence for BCinj1 inflow condition

        + default value is 1e-5

    - **'inj1_newton_nit'**: Newton iterations for BCinj1 inflow condition

        + default value is 10

    - **'psiroe'**: Harten correction parameter

        + default value is 0.1

    - **'prandtltb'**: turbulent Prandtl number (only active for 'model'='NSTurbulent')

        + default value is 0.92

    - **'ransmodel'**: RANS turbulence models. Possible values are

        + 'SA'      (Standard Spalart-Allmaras model)
        + 'SA_comp' (SA with mixing layer compressibility correction, see https://turbmodels.larc.nasa.gov/spalart.html)
        + default value is 'SA'
    
    - **'SA_add_RotCorr'**: add rotation correction to the SA model (SA-R)

        + 0
        + 1  (toggle rotation correction)
        + default value is 0

    - **'SA_add_LowRe'**: add low Reynolds correction to the SA model (SA-LRe)

        + 0
        + 1  (toggle low-Re correction)
        + default value is 0

    - **'DES'**: toggle the turbulence models for DES simulations. Possible values are

        + 'none'    (SA computation)
        + 'zdes1'   (ZDES mode 1, see https://link.springer.com/content/pdf/10.1007%2Fs00162-011-0240-z.pdf)
        + 'zdes1_w' (ZDES mode 1 by Chauvet)
        + 'zdes2'   (ZDES mode 2)
        + 'zdes2_w' (ZDES mode 2 by Chauvet)
        + 'zdes3'   (ZDES mode 3, see p118 https://tel.archives-ouvertes.fr/tel-01365361/document)
        + default value is 'none'

    - **'DES_debug'**: toggle data extraction for DES post-analysis. Possible values are

        + 0    
        + 1  (save delta and fd functions in the FlowSolution#Centers node)
        + default value is 'none'

    - **'sgsmodel'**: LES turbulence models. Possible values are

        + 'Miles'   (Monotonically integrated LES, for which ViscosityEddy==LaminarViscosity)
        + 'smsm'    (Selective mixed scale model, see Lenormand et al https://doi.org/10.1002/(SICI)1097-0363(20000229)32:4%3C369::AID-FLD943%3E3.0.CO;2-6)
        + default value is 'Miles'

    - **'extract_res'**: toggle residual extraction during simulations. Possible values are

        + 0  (don't save)
        + 1  (save div(F_Euler-F_viscous) in the FlowSolution#Centers node)
        + 2  (save dqdt + div(F_Euler-F_viscous) in the FlowSolution#Centers node)
        + default value is 0

    - **'source'**: toggle a source term for CFD simulations. Possible values are

        + 0   
        + 1  (read a source terme in the FlowSolution#Centers node. The conservative variables centers:Density_src, centers:MomentumX_src, centers:MomentumY_src, centers:MomentumZ_src and centers:EnergyStagnationDensity_src are used.)
        + default value is 0

    - **'ratiom'**: maximum cut-off for mut/mu

        + default value is 10000

    *Example of use:*

    * `Set numerics to zone (pyTree) <Examples/Fast/setNum2ZonesPT.py>`_:

    .. literalinclude:: ../../test/setNum2ZonesPT.py

----------------------------------------------------

.. py:function:: Fast.PyTree.loadTree(fileName='t.cgns', split='single', graph=False)

    Load a tree from a file (solution tree t, connectivity tree tc, or statistics tree ts). The file can be a single CGNS file (e.g. 't.cgns'), or can be split with one file per processor (e.g. 't/t_1.cgns', 't/t_2.cgns', etc.).
    
    When executed in sequential mode, the function returns a complete tree.

    When executed in mpi mode, the function returns a partial tree on each processor.

    .. warning:: in MPI mode, the tree must already be distributed with the correct number of processors (with 'proc' nodes within the zones).

    If graph is true, create and return the communication graph for Chimera and abutting transfers (for MPI use only). The dictionary is of the form of: graph={'graphID':[], 'graphIBC':[], 'procDict':[], 'procList':[]}.

    .. warning:: to create the communication graph, the loaded tree must be the connectivity tree (tc.cgns), not the solution tree (t.cgns).

    :param fileName: name of the file to load
    :type fileName: string
    :param split: file format: 'single' or 'multiple'
    :type split: string
    :param graph: true for returning communication graph
    :type graph: boolean
    :return: t or tc or ts or (tc, graph)

    * `Load single pyTree (pyTree) <Examples/Fast/loadTreePT.py>`_:

    .. literalinclude:: ../../test/loadTreePT.py

    .. note:: since Fast 4.1, loadFile and loadTree have been merged.

----------------------------------------------------

.. py:function:: Fast.PyTree.saveTree(t, fileName='restart.cgns', split='single', compress=0)
    
    Save a tree to a file (solution tree t, connectivity tree tc, or statistics tree ts). The file can be a single CGNS file (e.g. 't.cgns'), or can be split with one file per processor (e.g. 't/t_1.cgns', 't/t_2.cgns', etc.). The tree is also cleared of unnecessary flow fields and temporary work arrays.
    
    .. warning:: in MPI mode, the tree must already be distributed with the correct number of processors (with 'proc' nodes within the zones).

    Before saving the file, the flow fields can be compressed using the compress parameter.
    
    :param t: tree to save
    :type t: tree
    :param fileName: name of the file to load
    :type fileName: string
    :param split: file format: 'single' or 'multiple'
    :type split: string
    :param compress: compression parameter: 0 (none), 1 (cart only), or 2 (all)
    :type compress: integer
    
    * `Save single pyTree (pyTree) <Examples/Fast/saveTreePT.py>`_:

    .. literalinclude:: ../../test/saveTreePT.py

    .. note:: since Fast 4.1, saveFile and saveTree have been merged.

----------------------------------------------------

.. py:function:: Fast.PyTree.load(fileName='t.cgns', fileNameC=None, fileNameS=None, split='single')

    Load all trees at once (solution tree t, connectivity tree tc, and/or statistics tree ts), using multiple calls of the Fast.PyTree.loadTree function. The files can be single CGNS files (e.g. 't.cgns'), or can be split with one file per processor (e.g. 't/t_1.cgns', 't/t_2.cgns', etc.).

    .. warning:: in MPI mode, the trees must already be distributed with the correct number of processors (with 'proc' nodes within the zones).
    
    If fileNameC is not None, the connectivity tree is loaded and the communication graph is automatically created.

    If filenameS is not None, the statistics tree is loaded.

    :param fileName: name of the solution file to load
    :type fileName: string
    :param fileNameC: name of the connectivity file to load
    :type fileNameC: string
    :param fileNameS: name of the statistics file to load
    :type fileNameS: string
    :param split: file format: 'single' or 'multiple'
    :type split: string
    :return: (t, tc, ts, graph)

    * `Load pyTrees (pyTree) <Examples/Fast/loadPT.py>`_:

    .. literalinclude:: ../../test/loadPT.py


------------------------------------------------------------------

.. py:function:: Fast.PyTree.save(t, fileName='restart.cgns', tc=None, fileNameC='tc_restart.cgns', ts=None, fileNameS='tstat_restart.cgns', split='single', compress=0)
    
    Save all trees at once (solution tree t, connectivity tree tc, and/or statistics tree ts), using multiple calls of the Fast.PyTree.saveTree function. The files can be single CGNS files (e.g. 't.cgns'), or can be split with one file per processor (e.g. 't/t_1.cgns', 't/t_2.cgns', etc.). The trees are also cleared of unnecessary flow fields and temporary work arrays.

    .. warning:: in MPI mode, the trees must already be distributed with the correct number of processors (with 'proc' nodes within the zones).
    
    If tc is not None, the connectivity tree is saved.

    If ts is not None, the statistics tree is saved.

    Before saving the file, the flow fields can be compressed using the compress parameter.
    
    :param t: updated solution tree to save
    :type t: tree
    :param fileName: name of the solution file
    :type fileName: string
    :param tc: updated connectivity tree to save
    :type tc: tree
    :param fileNameC: name of the connectivity file
    :type fileNameC: string
    :param ts: updated statistics tree to save
    :type ts: tree
    :param fileNameS: name of the statistics file
    :type fileNameS: string
    :param split: file format: 'single' or 'multiple'
    :type split: string
    :param compress: compression parameter: 0 (none), 1 (cart only), or 2 (all)
    :type compress: integer

    *Example of use:*

    * `Save pyTrees (pyTree) <Examples/Fast/savePT.py>`_:

    .. literalinclude:: ../../test/savePT.py


.. toctree::
   :maxdepth: 2   


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

