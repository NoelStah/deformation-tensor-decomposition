# Deformation Tensor Tools
A Python package for analysing and decomposing deformation tensors using SymPy. This package includes utilities for computing the rate-of-strain tensor, rotation tensor, and more, with results easily displayed in LaTeX format.
***

## Repository Structure
```bash
.
├── README.md
├── deformation_tensor_tools
│   ├── __init__.py
│   └── deformation_tensor_tools.py
├── example_usage.ipynb
└── requirements.txt
```
## Features
- Compute the rate-of-strain tensor, its decomposition (volumetric and shear components) and its principal axes
- Compute the rotational tensor and the principal axis of rotation
- Perform a complete decomposition of the deformation tensor
- Easily display results in LaTeX using Jupyter Notebooks

## Installation
1. Clone the repository
```bash
git clone https://github.com/NoelStah/deformation-tensor-decomposition.git
```
2. Navigate to the project directory
```bash
cd deformation-tensor-decomposition
```
3. Install dependencies
This can be easily done with the help of the `requirements.txt` file
```bash
pip install -r requirements.txt
```
or manually by installing the SymPy package
```bash
pip install sympy
```

## Usage
1. Import the package
Ensure the package is accessible by including it in your python path. Import the SymPy package and the deformation tensor tools as follows:
```python
import sympy as sp
import deformation_tensor_tools as dtools
```

2. Example workflow
Refer to the example_usage.ipynb file for a detailed demonstration. Here is a basic example:
```python
import sympy as sp
import deformation_tensor_tools as dtools

# Define symbolic parameters
alpha, beta, gamma = sp.Symbols(r'\alpha \beta \gamma')

# Define a deformation tensor
deformation_tensor = sp.Matrix([
                                [alpha, beta, 0],
                                [0, 0, gamma],
                                [0, 0, 0]
                                ])

# Perform a complete decomposition
decomposition = dtools.complete_decomposition(deformation_tensor)

# Access components
rate_of_strain_tensor = decomposition[0]
rotation_tensor = decomposition[4]

# Display components
display(rate_of_strain_tensor)
display(rotation_tensor)
```
