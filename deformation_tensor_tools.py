# Helps calculating decompositions of the infinitesimal
# deformation tensor defined by del u_i/del x_j.
#######################################################

import sympy as sp

def calculate_dimension(deformation_tensor):
    '''
    If the deformation tensor is a square matrix, this function determines
    the dimension of it.
    '''
    # Raise error if the provided matrix is not a square matrix
    if (deformation_tensor.shape[0] != deformation_tensor.shape[1]):
        raise ValueError('Deformation tensor must be a square matrix.')
    else:
        dimension = deformation_tensor.shape[0]

    return(dimension)

def rate_of_strain_tensor(deformation_tensor):
    '''
    Calculates the rate of strain tensor 'epsilon' for a given deformation
    tensor.
    '''
    # Check for square matrix
    calculate_dimension(deformation_tensor)

    # Calculate the rate-of-strain tensor
    epsilon = (deformation_tensor + deformation_tensor.T) / sp.Rational(2) 

    return(epsilon)

def rate_of_strain_tensor_principle_axes(deformation_tensor):
    '''
    Calculates the principle axes of the rate of strain tensor 'epsilon'.
    '''
    # Check for square matrix
    calculate_dimension(deformation_tensor)

    # Calculate the rate of strain tensor
    epsilon = rate_of_strain_tensor(deformation_tensor)

    # Calculate the principle axis of rate-of-strain tensor
    epsilon_eigenvectors = epsilon.eigenvects()
    epsilon_axes = [eigvec for _, _, eigvecs in epsilon_eigenvectors
                    for eigvec in eigvecs]
    return(epsilon_axes)

def rate_of_strain_tensor_volumetric_part(deformation_tensor):
    '''
    Calculates the volumetric part of the rate of strain tensor 'epsilon'.
    '''
    # Check for square matrix
    dimension = calculate_dimension(deformation_tensor)

    # Calculate the rate of strain tensor
    epsilon = rate_of_strain_tensor(deformation_tensor)

    # Calculate the volumetric part of the rate-of-strain tensor
    epsilon_v = epsilon.trace() / sp.Rational(3) * sp.eye(dimension)

    return(epsilon_v)

def rate_of_strain_tensor_shear_part(deformation_tensor):
    '''
    Calculates the shear part of the rate of strain tensor 'epsilon'.
    '''
    # Check for square matrix
    calculate_dimension(deformation_tensor)

    # Calculate the rate of strain tensor and its volumetric part
    epsilon = rate_of_strain_tensor(deformation_tensor)
    epsilon_v = rate_of_strain_tensor_volumetric_part(deformation_tensor)

    # Calculate the shear part of the rate-of-strain tensor
    epsilon_s   = epsilon - epsilon_v

    return(epsilon_s)

def rotation_tensor(deformation_tensor):
    '''
    Calculates the rotation tensor 'omega' for a given deformation tensor.
    '''
    # Check for square matrix
    calculate_dimension(deformation_tensor)

    # Calculate the rotation tensor
    omega = (deformation_tensor - deformation_tensor.T) / sp.Rational(2)

    return(omega)

def principle_axis_of_rotation(deformation_tensor):
    '''
    Calculates the principle axis of rotation for a given deformation
    tensor (if it is a 3x3 matrix).
    '''
    # Check for square matrix
    dimension = calculate_dimension(deformation_tensor)

    # Calculate rotation tensor
    omega = rotation_tensor(deformation_tensor)

    # Calculate the principle axis of rotation (for 3x3 matrices)
    if (dimension != 3):
        raise ValueError('Deformation tensor must be a 3x3 matrix.')
    else:
        rotational_axis = sp.Matrix([
                                    [(omega[1,2]-omega[2,1]) / sp.Rational(2)],
                                    [(omega[2,0]-omega[0,2]) / sp.Rational(2)],
                                    [(omega[0,1]-omega[1,0]) / sp.Rational(2)]
                                    ])

    return(rotational_axis)

def complete_decomposition(deformation_tensor):
    '''
    Completely decomposes a given deformation tensor into its:
        - rate-of-strain tensor
        - principle axes of the rate-of-strain tensor
        - volumetric part of the rate-of-strian tensor
        - shear part of the rate-of-strain tensor
        - rotation tensor
        - principle axis of rotation (for a 3x3 matrix)
    '''
    # Check for square matrix
    dimension = calculate_dimension(deformation_tensor)

    # Calculate the rate-of-strain tensor
    epsilon = rate_of_strain_tensor(deformation_tensor)

    # Calculate the principle axis of rate-of-strain tensor
    epsilon_axes = rate_of_strain_tensor_principle_axes(deformation_tensor)

    # Calculate the volumetric part of the rate-of-strain tensor
    epsilon_v = rate_of_strain_tensor_volumetric_part(deformation_tensor)

    # Calculate the shear part of the rate-of-strain tensor
    epsilon_s = rate_of_strain_tensor_shear_part(deformation_tensor)

    # Calculate the rotation tensor
    omega = rotation_tensor(deformation_tensor)

    # Calculate the principle axis of rotation (for 3x3 matrices)
    rotational_axis = principle_axis_of_rotation(deformation_tensor)
    return(epsilon, epsilon_axes, epsilon_v, epsilon_s, omega, rotational_axis)
