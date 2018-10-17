.. Fast documentation master file

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

   Fast.PyTree.setNum2Base
   Fast.PyTree.setNum2Zones
   Fast.PyTree.load
   Fast.PyTree.save
   Fast.PyTree.loadFile
   Fast.PyTree.saveFile
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

    - **'temporal_scheme'**: possible values are

        + 'explicit' (RK3 scheme, see p49 http://publications.onera.fr/exl-php/util/documents/accede_document.php )
        + 'implicit' (BDF2 or BDF1 if local time stepping  )
        + 'implicit_local' (see p107 http://publications.onera.fr/exl-php/docs/ILS_DOC/227155/DOC356618_s1.pdf)
        + default value is 'implicit'

    - **'ss_iteration'**: 

        + Newton Iterations for implicit schemes
        + default value is 30

    - **'modulo_verif'**: 

        + period of computation for: cfl (RK3 or BDF2), newton convergence (all temporal_scheme) and predictor estimation for 'implicit_local' scheme
        + default value is 200

    *Example of use:*

    * `Set numerics to base (pyTree) <Examples/Fast/setNum2BasePT.py>`_:

    .. literalinclude:: ../build/Examples/Fast/setNum2BasePT.py

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

    - **'scheme'**: possible values are
 
        + 'ausmpred' (see  p49 https://tel.archives-ouvertes.fr/pastel-00834850/document)
        + 'sensor'   (see  p50 https://tel.archives-ouvertes.fr/pastel-00834850/document)
        + 'roe'
        + default value is 'ausmpred'

    - **'slope'**: possible values are

        + 'o3'      (third order, see  p50 https://tel.archives-ouvertes.fr/pastel-00834850/document)
        + 'o1'      (first order,  only valid for roe scheme)
        + 'minmod'  (only valid for roe scheme)
        + default value is 'o3'

    - **'motion'**: possible values are

        + 'none' (no motion)
        + 'rigid' (ALE without deformation see p47 https://tel.archives-ouvertes.fr/tel-01011273/document)
        + default value is 'none'

    - **'time_step'**: 

        + value of time step
        + default value is 1e-4

    - **'time_step_nature'**: 

        + 'global'
        + 'local'
        + default value is 'global'

    - **'epsi_newton'**: 

        + newton stopping criteria on Loo norm
        + default value is 0.1
    
    - **'inj1_newton_tol'**:

        + Newton tolerence for BCinj1 inflow condition
        + default value is 1e-5

    - **'inj1_newton_nit'**: 

        + Newton Iteration for BCinj1 inflow condition
        + default value is 10

    - **'psiroe'**: 

        + Harten correction
        + default value is 0.1
        
    - **'cfl'**: 

        + usefull only if 'time_step_nature'='local'
        + default value is 1

    - **'model'**: possible values are
 
        + 'Euler'
        + 'NSLaminar'
        + 'NSTurbulent'(only Spalart available)
        + default value is 'Euler'

    - **'prandtltb'**:

        + turbulent Prandtl number (only active for 'model'='NSTurbulent')
        + default value is 0.92

    - **'ransmodel'**: possible values are

        + "SA"      (Standard Spalart-Allmaras model)
        + "SA_comp" (SA with mixing layer compressible correction https://turbmodels.larc.nasa.gov/spalart.html)
        + default value is 'SA'

    - **'DES'**: possible values are

        + "none"    (SA computation)
        + "zdes1"   (mode 1, https://link.springer.com/content/pdf/10.1007%2Fs00162-011-0240-z.pdf)
        + "zdes1_w" (mode 1 by Chauvet)
        + "zdes2"   (mode 2)
        + "zdes2_w" (mode 2 by Chauvet)
        + "zdes3"   (mode 3, see p118 https://tel.archives-ouvertes.fr/tel-01365361/document))
        + default value is 'none'

    - **'DES_debug'**: possible values are

        + "none"    
        + "active"  (save delta and fd functions in the FlowSolution#Centers node)
        + default value is 'none'

    - **'sgsmodel'**: possible values are

        + "Miles"   (ViscosityEddy==LaminarViscosity)
        + "smsm"    (Selective Mixed Scale model, Lenormand et al, (2000), LES of sub and supersonic channel flow at moderate Re. Int. J. Numer. Meth. Fluids, 32: 369–406)
        + default value is 'Miles'

    - **'extract_res'**: possible values are

        + 0   
        + 1  (save div(F_Euler-F_viscous) in the FlowSolution#Centers node)
        + default value is 0

    *Example of use:*

    * `Set numerics to zone (pyTree) <Examples/Fast/setNum2ZonesPT.py>`_:

    .. literalinclude:: ../build/Examples/Fast/setNum2ZonesPT.py

..  ----------------------------------------------------
.. ..py:function:: Fast.PyTree.metric(a)
.. Compute the metric needed by the solvers.
.. :param a: input data
.. :type a: Zone, Base, pyTree
.. :return: a, metrics
.. * `Compute metric for zones (pyTree) <Examples/Fast/metricPT.py>`_:
.. .. literalinclude:: ../build/Examples/Fast/metricPT.py


----------------------------------------------------

.. py:function:: Fast.PyTree.load(fileName='t.cgns', fileNameC= 'tc.cgns', fileNameS='tstat.cgns', split='single', NP=0)
    
    Load computation tree t from file. 
    Optionaly load tc (connectivity file) or tstat (statistics file).
    Returns also the graph as a dictionary {'graphID', 'graphIBC', 'procDict', 'procList'}.
    If you run in mpi, NP must be the number of processor.
    If you run in seq mode, NP must be 0 or a negative number.
    If split='single', a single file is written.
    If split='multiple', different files are created
    depending on the proc number of each zone (restart/restart_0.cgns, ...).

    :param a: input data
    :type a: pyTree
    :param fileName: name of file for save
    :type fileName: string
    :param split: 'single' or 'multiple'
    :type split: string
    :param NP: number of processors
    :type NP: int
    :return: t, tc, ts, graph
    :rtype: tuple

    * `Load pyTree (pyTree) <Examples/Fast/loadPT.py>`_:

    .. literalinclude:: ../build/Examples/Fast/loadPT.py


------------------------------------------------------------------

.. py:function:: Fast.PyTree.save(t, fileName='restart.cgns', split='single', temporal_scheme='implicit', NP=0)
    
    Save computation tree t in file. 
    If you run in mpi, NP must be the number of processor.
    If you run in seq mode, NP must be 0 or a negative number.
    If split='single', a single file is written.
    If split='multiple', different files are created
    depending on the proc number of each zone (restart/restart_0.cgns, ...).
    
    :param a: input data
    :type a: pyTree
    :param fileName: name of file for save
    :type fileName: string
    :param split: 'single' or 'multiple'
    :type split: string
    :param NP: number of processors
    :type NP: int

    *Example of use:*

    * `Save pyTree (pyTree) <Examples/Fast/savePT.py>`_:

    .. literalinclude:: ../build/Examples/Fast/savePT.py

----------------------------------------------------

.. py:function:: Fast.PyTree.loadFile(fileName='t.cgns', split='single', mpirun=False)
    
    Load tree from file. The tree must be already distributed (with 'proc' nodes).
    The file can be a single CGNS file ("t.cgns") or a splitted per processor 
    CGNS file ("t/t_1.cgns", "t/t_2.cgns", ...)
    
    If you run in sequential mode, mpirun must be false.
    The function returns a full tree.

    If you run in mpi mode, mpirun must be true.
    The function returns a partial tree on each processor.

    :param fileName: name of file for load
    :type fileName: string
    :param split: 'single' or 'multiple'
    :type split: string
    :param mpirun: true if python is run with mpirun
    :type mpirun boolean
    :return: t
    :rtype: CGNS tree

    * `Load single pyTree (pyTree) <Examples/Fast/loadFilePT.py>`_:

    .. literalinclude:: ../build/Examples/Fast/loadFilePT.py

----------------------------------------------------

.. py:function:: Fast.PyTree.saveFile(fileName='t.cgns', split='single', mpirun=False)
    
    Save tree to file. The tree must be already distributed (with 'proc' nodes).

    The file can be a single CGNS file ("t.cgns") or a splitted per processor 
    CGNS file ("t/t_1.cgns", "t_2.cgns", ...)
    
    If you run in seq mode, mpirun must be false.
    
    If you run in mpi mode, mpirun must be true.
    
    :param fileName: name of file for load
    :type fileName: string
    :param split: 'single' or 'multiple'
    :type split: string
    :param mpirun: true if python is run with mpirun
    :type mpirune: boolean
    
    * `Save single pyTree (pyTree) <Examples/Fast/saveFilePT.py>`_:

    .. literalinclude:: ../build/Examples/Fast/saveFilePT.py


.. toctree::
   :maxdepth: 2   


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

