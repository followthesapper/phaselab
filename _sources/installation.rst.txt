Installation
============

PhaseLab can be installed via pip with various optional dependencies.

Basic Installation
------------------

.. code-block:: bash

   pip install phaselab

This installs the core package with essential dependencies.

Optional Dependencies
---------------------

**With Quantum Support (IBM Quantum)**

.. code-block:: bash

   pip install phaselab[quantum]

Includes Qiskit and IBM Quantum Runtime for hardware validation.

**With ATLAS-Q Integration**

.. code-block:: bash

   pip install phaselab[atlas]

Enables ATLAS-Q tensor network backend for accelerated coherence computation.

**With Plotting Support**

.. code-block:: bash

   pip install phaselab[plotting]

Includes matplotlib and seaborn for visualization.

**Full Installation**

.. code-block:: bash

   pip install phaselab[all]

Installs all optional dependencies.

Development Installation
------------------------

For contributing or development:

.. code-block:: bash

   git clone https://github.com/followthesapper/phaselab.git
   cd phaselab
   pip install -e ".[dev]"

This installs:

- All optional dependencies
- Testing tools (pytest, coverage)
- Documentation tools (sphinx)
- Code quality tools (black, ruff, mypy)

System Requirements
-------------------

- Python 3.9 or higher
- NumPy 1.21+
- SciPy 1.7+
- Pandas 1.3+

**For Quantum Features:**

- Qiskit 1.0+ (for IBM Quantum)
- qiskit-ibm-runtime (for hardware access)
- IBM Quantum account (free tier available)

**For ATLAS-Q Acceleration:**

- PyTorch 2.0+ (optional, for GPU)
- CUDA 11.8+ (optional, for GPU acceleration)

Verifying Installation
----------------------

.. code-block:: python

   import phaselab as pl

   # Check version
   print(f"PhaseLab version: {pl.__version__}")

   # Test core functionality
   R_bar = pl.coherence_score([0.1, 0.15, 0.12])
   print(f"Coherence: {R_bar:.4f}")
   print(f"Status: {pl.go_no_go(R_bar)}")

   # Check ATLAS-Q availability
   from phaselab.quantum import is_atlas_q_available
   print(f"ATLAS-Q available: {is_atlas_q_available()}")

Expected output:

.. code-block:: text

   PhaseLab version: 0.6.1
   Coherence: 0.9995
   Status: GO
   ATLAS-Q available: True/False

Troubleshooting
---------------

**ImportError: No module named 'phaselab'**

Ensure you've installed the package in the correct Python environment:

.. code-block:: bash

   python -m pip install phaselab
   python -c "import phaselab; print(phaselab.__version__)"

**ATLAS-Q not available**

ATLAS-Q requires separate installation:

.. code-block:: bash

   pip install atlas-quantum

Or install from source:

.. code-block:: bash

   git clone https://github.com/followthesapper/ATLAS-Q.git
   cd ATLAS-Q
   pip install -e .

**IBM Quantum authentication errors**

Set your IBM Quantum token:

.. code-block:: python

   import os
   os.environ['IBM_QUANTUM_TOKEN'] = 'your-token-here'

Or save credentials:

.. code-block:: python

   from qiskit_ibm_runtime import QiskitRuntimeService
   QiskitRuntimeService.save_account(channel="ibm_quantum", token="your-token")
