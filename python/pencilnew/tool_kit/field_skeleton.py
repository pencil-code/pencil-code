# field_skeleton.py
# Written by Simon Candelaresi (iomsn1@gmail.com)

"""
Finds the structure of the field's skeleton, i.e. null points, separatrix
layers and separators using the trilinear method
Haynes-Parnell-2007-14-8-PhysPlasm (http://dx.doi.org/10.1063/1.2756751).
"""

import numpy as np


class NullPoint(object):
    """
    Contains the position of the null points.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.null = []

    def find_nullpoints(self, var, field):
        """
        Find the null points to the field 'field' with information from 'var'.
        """

        def __triLinear_interpolation(x, y, z, coefTri):
            """ Compute the interpolated field at (normalized) x, y, z."""
            return coefTri[0] + coefTri[1]*x + coefTri[2]*y + coefTri[3]*x*y +\
                   coefTri[4]*z + coefTri[5]*x*z + coefTri[6]*y*z + \
                   coefTri[7]*x*y*z

        def __gradB1(x, y, z, coefTri, dd):
            """ Compute the inverse of the gradient of B. """
            gb1 = np.zeros((3, 3))
            gb1[0, :] = (__triLinear_interpolation(x+dd, y, z, coefTri) - \
                         __triLinear_interpolation(x-dd, y, z, coefTri))/(2*dd)
            gb1[1, :] = (__triLinear_interpolation(x, y+dd, z, coefTri) - \
                         __triLinear_interpolation(x, y-dd, z, coefTri))/(2*dd)
            gb1[2, :] = (__triLinear_interpolation(x, y, z+dd, coefTri) - \
                         __triLinear_interpolation(x, y, z-dd, coefTri))/(2*dd)
            # Invert the matrix.
            if np.linalg.det(gb1) != 0 and not np.max(np.isnan(gb1)):
                gb1 = np.matrix(gb1).I
            else:
                gb1 *= 0
            return gb1

        def __newton_raphson(xyz0, coefTri):
            """ Newton-Raphson method for finding null-points. """
            xyz = np.array(xyz0)
            iterMax = 10
            dd = 1e-4
            tol = 1e-5
            
            for i in range(iterMax):
                diff = __triLinear_interpolation(xyz[0], xyz[1],
                                                 xyz[2], coefTri) * \
                       __gradB1(xyz[0], xyz[1], xyz[2], coefTri, dd)
                diff = np.array(diff)[0]
                xyz = xyz - diff
                if any(abs(diff) < tol) or any(abs(diff) > 1):
                    return xyz

            return np.array(xyz)[0]

        # 1) Reduction step.
        # Find all cells for which all three field components change sign.
        sign_field = np.sign(field)
        reduced_cells = True
        for comp in range(3):
            reduced_cells *= \
            (sign_field[comp, 1:, 1:, 1:]*sign_field[comp, :-1, 1:, 1:] < 0) + \
            (sign_field[comp, 1:, 1:, 1:]*sign_field[comp, 1:, :-1, 1:] < 0) + \
            (sign_field[comp, 1:, 1:, 1:]*sign_field[comp, 1:, 1:, :-1] < 0) + \
            (sign_field[comp, 1:, 1:, 1:]*sign_field[comp, :-1, :-1, 1:] < 0) + \
            (sign_field[comp, 1:, 1:, 1:]*sign_field[comp, 1:, :-1, :-1] < 0) + \
            (sign_field[comp, 1:, 1:, 1:]*sign_field[comp, :-1, 1:, :-1] < 0) + \
            (sign_field[comp, 1:, 1:, 1:]*sign_field[comp, :-1, :-1, :-1] < 0)

        # Find null points in these cells.
        for cell_idx in range(np.sum(reduced_cells)):
            # 2) Analysis step.

            # Find the indices of the cell where to look for the null point.
            idx_x = np.where(reduced_cells == True)[2][cell_idx]
            idx_y = np.where(reduced_cells == True)[1][cell_idx]
            idx_z = np.where(reduced_cells == True)[0][cell_idx]
            x = var.x
            y = var.y
            z = var.z

            # Compute the coefficients for the trilinear interpolation.
            coefTri = np.zeros((8, 3))
            coefTri[0] = field[:, idx_z, idx_y, idx_x]
            coefTri[1] = field[:, idx_z, idx_y, idx_x+1] - \
                         field[:, idx_z, idx_y, idx_x]
            coefTri[2] = field[:, idx_z, idx_y+1, idx_x] - \
                         field[:, idx_z, idx_y, idx_x]
            coefTri[3] = field[:, idx_z, idx_y+1, idx_x+1] - \
                         field[:, idx_z, idx_y, idx_x+1] - \
                         field[:, idx_z, idx_y+1, idx_x] + \
                         field[:, idx_z, idx_y, idx_x]
            coefTri[4] = field[:, idx_z+1, idx_y, idx_x] - \
                         field[:, idx_z, idx_y, idx_x]
            coefTri[5] = field[:, idx_z+1, idx_y, idx_x+1] - \
                         field[:, idx_z, idx_y, idx_x+1] - \
                         field[:, idx_z+1, idx_y, idx_x] + \
                         field[:, idx_z, idx_y, idx_x]
            coefTri[6] = field[:, idx_z+1, idx_y+1, idx_x] - \
                         field[:, idx_z, idx_y+1, idx_x] - \
                         field[:, idx_z+1, idx_y, idx_x] + \
                         field[:, idx_z, idx_y, idx_x]
            coefTri[7] = field[:, idx_z+1, idx_y+1, idx_x+1] - \
                         field[:, idx_z, idx_y+1, idx_x+1] - \
                         field[:, idx_z+1, idx_y, idx_x+1] - \
                         field[:, idx_z+1, idx_y+1, idx_x] + \
                         field[:, idx_z, idx_y, idx_x+1] + \
                         field[:, idx_z, idx_y+1, idx_x] + \
                         field[:, idx_z+1, idx_y, idx_x] - \
                         field[:, idx_z, idx_y, idx_x]

            # Find the intersection of the curves field_i = field_j = 0.
            # The units are first normalized to the unit cube from 000 to 111.

            # face 1
            intersection = False
            coefBi = np.zeros((4, 3))
            coefBi[0] = coefTri[0] + coefTri[4]*0
            coefBi[1] = coefTri[1] + coefTri[5]*0
            coefBi[2] = coefTri[2] + coefTri[6]*0
            coefBi[3] = coefTri[3] + coefTri[7]*0
            # Find the roots for x0 and y0.
            polynomial = np.array([coefBi[1, 0]*coefBi[3, 1] -
                                   coefBi[1, 1]*coefBi[3, 0],
                                   coefBi[0, 0]*coefBi[3, 1] +
                                   coefBi[1, 0]*coefBi[2, 1] -
                                   coefBi[0, 1]*coefBi[3, 0] -
                                   coefBi[2, 0]*coefBi[1, 1],
                                   coefBi[0, 0]*coefBi[2, 1] -
                                   coefBi[0, 1]*coefBi[2, 0]])
            roots_x = np.roots(polynomial)
            if len(roots_x) == 0:
                roots_x = -np.ones(2)
            if ((roots_x[0] >= 0) and (roots_x[0] <= 1)) or \
            ((roots_x[1] >= 0) and (roots_x[1] <= 1)):
                roots_y = -(coefBi[0, 0] + coefBi[1, 0]*roots_x)/ \
                           (coefBi[2, 0] + coefBi[3, 0]*roots_x)
                if ((roots_y[0] >= 0) and (roots_y[0] <= 1)) or \
                ((roots_y[1] >= 0) and (roots_y[1] <= 1)):
                    intersection = True
            if intersection:
                xyz0 = [roots_x[1], roots_y[1], 0]
                xyz = __newton_raphson(xyz0, coefTri)
                xyz1 = [xyz[0]*var.dx + x[idx_x],
                       xyz[1]*var.dy + y[idx_y],
                       xyz[2]*var.dz + z[idx_z]]
                print xyz1

            # face 2
            intersection = False
            coefBi = np.zeros((4, 3))
            coefBi[0] = coefTri[0] + coefTri[4]*1
            coefBi[1] = coefTri[1] + coefTri[5]*1
            coefBi[2] = coefTri[2] + coefTri[6]*1
            coefBi[3] = coefTri[3] + coefTri[7]*1
            # Find the roots for x0 and y0.
            polynomial = np.array([coefBi[1, 0]*coefBi[3, 1] -
                                   coefBi[1, 1]*coefBi[3, 0],
                                   coefBi[0, 0]*coefBi[3, 1] +
                                   coefBi[1, 0]*coefBi[2, 1] -
                                   coefBi[0, 1]*coefBi[3, 0] -
                                   coefBi[2, 0]*coefBi[1, 1],
                                   coefBi[0, 0]*coefBi[2, 1] -
                                   coefBi[0, 1]*coefBi[2, 0]])
            roots_x = np.roots(polynomial)
            if len(roots_x) == 0:
                roots_x = -np.ones(2)
            if ((roots_x[0] >= 0) and (roots_x[0] <= 1)) or \
            ((roots_x[1] >= 0) and (roots_x[1] <= 1)):
                roots_y = -(coefBi[0, 0] + coefBi[1, 0]*roots_x)/ \
                           (coefBi[2, 0] + coefBi[3, 0]*roots_x)
                if ((roots_y[0] >= 0) and (roots_y[0] <= 1)) or \
                ((roots_y[1] >= 0) and (roots_y[1] <= 1)):
                    intersection = True
            if intersection:
                xyz0 = [roots_x[1], roots_y[1], 0]
                xyz = __newton_raphson(xyz0, coefTri)
                xyz1 = [xyz[0]*var.dx + x[idx_x],
                       xyz[1]*var.dy + y[idx_y],
                       xyz[2]*var.dz + z[idx_z]]
                print xyz1

            # face 3
            intersection = False
            coefBi = np.zeros((4, 3))
            coefBi[0] = coefTri[0] + coefTri[2]*0
            coefBi[1] = coefTri[1] + coefTri[3]*0
            coefBi[2] = coefTri[4] + coefTri[6]*0
            coefBi[3] = coefTri[5] + coefTri[7]*0
            # Find the roots for x0 and z0.
            polynomial = np.array([coefBi[1, 0]*coefBi[3, 2] -
                                   coefBi[1, 2]*coefBi[3, 0],
                                   coefBi[0, 0]*coefBi[3, 2] +
                                   coefBi[1, 0]*coefBi[2, 2] -
                                   coefBi[0, 2]*coefBi[3, 0] -
                                   coefBi[2, 0]*coefBi[1, 2],
                                   coefBi[0, 0]*coefBi[2, 2] -
                                   coefBi[0, 2]*coefBi[2, 0]])
            roots_x = np.roots(polynomial)
            if len(roots_x) == 0:
                roots_x = -np.ones(2)
            if ((roots_x[0] >= 0) and (roots_x[0] <= 1)) or \
            ((roots_x[1] >= 0) and (roots_x[1] <= 1)):
                roots_z = -(coefBi[0, 0] + coefBi[1, 0]*roots_x)/ \
                           (coefBi[2, 0] + coefBi[3, 0]*roots_x)
                if ((roots_z[0] >= 0) and (roots_z[0] <= 1)) or \
                ((roots_z[1] >= 0) and (roots_z[1] <= 1)):
                    intersection = True
            if intersection:
                xyz0 = [roots_x[1], 0, roots_z[1]]
                xyz = __newton_raphson(xyz0, coefTri)
                xyz1 = [xyz[0]*var.dx + x[idx_x],
                       xyz[1]*var.dy + y[idx_y],
                       xyz[2]*var.dz + z[idx_z]]
                print xyz1

            # face 4
            intersection = False
            coefBi = np.zeros((4, 3))
            coefBi[0] = coefTri[0] + coefTri[2]*1
            coefBi[1] = coefTri[1] + coefTri[3]*1
            coefBi[2] = coefTri[4] + coefTri[6]*1
            coefBi[3] = coefTri[5] + coefTri[7]*1
            # Find the roots for x0 and z0.
            polynomial = np.array([coefBi[1, 0]*coefBi[3, 2] -
                                   coefBi[1, 2]*coefBi[3, 0],
                                   coefBi[0, 0]*coefBi[3, 2] +
                                   coefBi[1, 0]*coefBi[2, 2] -
                                   coefBi[0, 2]*coefBi[3, 0] -
                                   coefBi[2, 0]*coefBi[1, 2],
                                   coefBi[0, 0]*coefBi[2, 2] -
                                   coefBi[0, 2]*coefBi[2, 0]])
            roots_x = np.roots(polynomial)
            if len(roots_x) == 0:
                roots_x = -np.ones(2)
            if ((roots_x[0] >= 0) and (roots_x[0] <= 1)) or \
            ((roots_x[1] >= 0) and (roots_x[1] <= 1)):
                roots_z = -(coefBi[0, 0] + coefBi[1, 0]*roots_x)/ \
                           (coefBi[2, 0] + coefBi[3, 0]*roots_x)
                if ((roots_z[0] >= 0) and (roots_z[0] <= 1)) or \
                ((roots_z[1] >= 0) and (roots_z[1] <= 1)):
                    intersection = True
            if intersection:
                xyz0 = [roots_x[1], 0, roots_z[1]]
                xyz = __newton_raphson(xyz0, coefTri)
                print "shape = ", x.shape, xyz.shape
#                xyz1 = [xyz[0]*var.dx + x[idx_x],
#                       xyz[1]*var.dy + y[idx_y],
#                       xyz[2]*var.dz + z[idx_z]]
#                print xyz1

            # face 5
            intersection = False
            coefBi = np.zeros((4, 3))
            coefBi[0] = coefTri[0] + coefTri[1]*0
            coefBi[1] = coefTri[2] + coefTri[3]*0
            coefBi[2] = coefTri[4] + coefTri[5]*0
            coefBi[3] = coefTri[6] + coefTri[7]*0
            # Find the roots for y0 and z0.
            polynomial = np.array([coefBi[1, 1]*coefBi[3, 2] -
                                   coefBi[1, 2]*coefBi[3, 1],
                                   coefBi[0, 1]*coefBi[3, 2] +
                                   coefBi[1, 1]*coefBi[2, 2] -
                                   coefBi[0, 2]*coefBi[3, 1] -
                                   coefBi[2, 1]*coefBi[1, 2],
                                   coefBi[0, 1]*coefBi[2, 2] -
                                   coefBi[0, 2]*coefBi[2, 1]])
            roots_y = np.roots(polynomial)
            if len(roots_y) == 0:
                roots_y = -np.ones(2)
            if ((roots_y[0] >= 0) and (roots_y[0] <= 1)) or \
            ((roots_y[1] >= 0) and (roots_y[1] <= 1)):
                roots_z = -(coefBi[0, 1]+coefBi[1, 1]*roots_x)/ \
                           (coefBi[2, 1]+coefBi[3, 1]*roots_x)
                if ((roots_z[0] >= 0) and (roots_z[0] <= 1)) or \
                ((roots_z[1] >= 0) and (roots_z[1] <= 1)):
                    intersection = True
            if intersection:
                xyz0 = [0, roots_y[1], roots_z[1]]
                xyz = __newton_raphson(xyz0, coefTri)
                xyz1 = [xyz[0]*var.dx + x[idx_x],
                       xyz[1]*var.dy + y[idx_y],
                       xyz[2]*var.dz + z[idx_z]]
                print xyz1

            # Convert from unit cube into real values.
            self.null = [xyz[0]*var.dx + x[idx_x],
                         xyz[1]*var.dy + y[idx_y],
                         xyz[2]*var.dz + z[idx_z]]
                         
            # face 6
            intersection = False
            coefBi = np.zeros((4, 3))
            coefBi[0] = coefTri[0] + coefTri[1]*1
            coefBi[1] = coefTri[2] + coefTri[3]*1
            coefBi[2] = coefTri[4] + coefTri[5]*1
            coefBi[3] = coefTri[6] + coefTri[7]*1
            # Find the roots for y0 and z0.
            polynomial = np.array([coefBi[1, 1]*coefBi[3, 2] -
                                   coefBi[1, 2]*coefBi[3, 1],
                                   coefBi[0, 1]*coefBi[3, 2] +
                                   coefBi[1, 1]*coefBi[2, 2] -
                                   coefBi[0, 2]*coefBi[3, 1] -
                                   coefBi[2, 1]*coefBi[1, 2],
                                   coefBi[0, 1]*coefBi[2, 2] -
                                   coefBi[0, 2]*coefBi[2, 1]])
            roots_y = np.roots(polynomial)
            if len(roots_y) == 0:
                roots_y = -np.ones(2)
            if ((roots_y[0] >= 0) and (roots_y[0] <= 1)) or \
            ((roots_y[1] >= 0) and (roots_y[1] <= 1)):
                roots_z = -(coefBi[0, 1]+coefBi[1, 1]*roots_x)/ \
                           (coefBi[2, 1]+coefBi[3, 1]*roots_x)
                if ((roots_z[0] >= 0) and (roots_z[0] <= 1)) or \
                ((roots_z[1] >= 0) and (roots_z[1] <= 1)):
                    intersection = True
            if intersection:
                xyz0 = [0, roots_y[1], roots_z[1]]
                xyz = __newton_raphson(xyz0, coefTri)
                xyz1 = [xyz[0]*var.dx + x[idx_x],
                       xyz[1]*var.dy + y[idx_y],
                       xyz[2]*var.dz + z[idx_z]]
                print xyz1

