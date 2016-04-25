# -*- coding: utf-8 -*-

# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2016 Alberto Fernandez                                  *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU Lesser General Public License (LGPL)    *
# *   as published by the Free Software Foundation; either version 2 of     *
# *   the License, or (at your option) any later version.                   *
# *   for detail see the LICENCE text file.                                 *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU Library General Public License for more details.                  *
# *                                                                         *
# *   You should have received a copy of the GNU Library General Public     *
# *   License along with this program; if not, write to the Free Software   *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
# *   USA                                                                   *
# ***************************************************************************

__author__ = "Alberto Fernandez"

import numpy as np
import xml.etree.ElementTree as ET

class Lamina:
    """ This class describes a lamina. A lamina has the properties defined on its on axis.
    Once it is laid up we use a layer class.
    This class keep information about the kind of lamina used. When properties are calculated, it is
    equivalent to using Puck's micromechanical model
    """

    def __init__(self, name='Name of lamina', uuid=''):
        """Initialization function creates a xml tree element to keep lamina properties. It follows a structure similar
        to the one seen in ElamX. In the future it could be a more independent structure"""
        self.xml = ET.Element('lamina')
        self.xml.set('name', name)
        self.xml.set('uuid', uuid)
        properties = ['E11', 'E22', 'G12', 'poisson12', 'poisson21', 't', 'rho', 'sigma1_t', 'sigma1_c',
                      'sigma2_t', 'sigma2_c', 'tau12']
        for prop in properties:
            self.xml.append(ET.Element(prop))
        self.Q = np.zeros((3, 3))

    def set_micro_prop(self, mfraction=0.5, matrix='epoxy', fiber='hs carbon', arealweight=300, arch='uni', angles=[0]):
        """This function takes care of a special initialization of the class in case it's properties are derived from
        the properties of its constituting materials. Most of the work is assigning either default values to the lamina
        or selecting appropriated ones from selected materials.
        Arch is the architecture of the lamina. It can be 'uni' for either unidirectional or multiaxial laminas.
        Multiaxial laminas have more than one angle in angles. For cloth use either 'plain', 'twill' or 'satin'"""
        # TODO there is a regresion on the way multiaxial reinforcements are treated. It shoul be refactored
        self.fiber = fiber.lower()
        self.matrix = matrix.lower()
        self.mfraction = mfraction
        self.arealweight = arealweight
        self.arch = arch.lower()
        self.angles = angles
        self.rho_f = 1740  # Default value for hs_carbon
        self.rho_m = 1200  # Default value for matrix
        self.Ex_f = 230000
        self.Ey_f = 28000
        self.Gxy_f = 50000
        self.poisson_f = 0.23
        self.E_m = 3600
        self.poisson_m = 0.35
        self.G_m = 0.0
        # Change values depending on materials. Densities in kg/m3
        if fiber == 'eglass':
            self.rho_f = 2540
            self.Ex_f = 73000
            self.Ey_f = 73000
            self.Gxy_f = 30000
            self.poisson_f = 0.18
        elif fiber == 'aramid':
            self.rho_f = 1440
            self.Ex_f = 124000
            self.Ey_f = 6900
            self.Gxy_f = 2800
            self.poisson_f = 0.36
        elif fiber == 'hm carbon':
            self.rho_f = 1810
            self.Ex_f = 392000
            self.Ey_f = 15000
            self.Gxy_f = 28600
            self.poisson_f = 0.20
        if matrix == 'polyester':
            self.rho_m = 1200
            self.E_m = 3000
            self.poisson_m = 0.316
        # calling the function to calculate properties from micromechanics and putting the list
        # returned into the xml tree
        properties = self.calc_local_properties()
        for prop in properties:
            self.xml.find(prop).text = str(properties[prop])
        self.xml.find('poisson21').text = str(self.calc_poisson21())
        return

    def calc_local_properties(self):
        """This function calculates the properties of a lamina from mass fraction and materials
        It uses default values for high strength carbon in an epoxy matrix. A mass fraction of 0.5
        is also assumed if the function is called without parameters"""
        self.vfraction = self.mfraction / (self.mfraction + (1 - self.mfraction) * self.rho_f / self.rho_m)
        E11 = self.vfraction * self.Ex_f + (1 - self.vfraction) * self.E_m
        E22 = self.E_m / (1 - self.poisson_m**2) * (1 + 0.85 * self.vfraction**2) / (
            (1 - self.vfraction)**1.25 + self.vfraction * self.E_m / (self.Ey_f * (1 - self.poisson_m**2)))
        poisson12 = self.vfraction * self.poisson_f + (1 - self.vfraction) * self.poisson_m
        self.G_m = self.E_m / (2 * (1 + self.poisson_m))
        G12 = self.G_m * ((1 + 0.8 * self.vfraction**0.8) / ((1 - self.vfraction)**1.25 + self.G_m /
                                                                  self.Gxy_f * self.vfraction))
        # units for thickness are mm. Input units for areal weight are g/m3
        #TODO cambiar estimacion de espesor a la clase ply
        self.t = self.arealweight * (1/self.rho_f + (1-self.mfraction)/(self.mfraction * self.rho_m))
        return {'E11': E11, 'E22': E22, 'G12': G12, 'poisson12': poisson12}

    def calc_poisson21(self):
        """Utility function to help with the calculations from outside the class"""
        E11 = float(self.xml.find('E11').text)
        E22 = float(self.xml.find('E22').text)
        poisson12 = float(self.xml.find('poisson12').text)
        return poisson12 * E22 / E11

    def calc_micro_Q_matrix(self):
        """Calculate matrix Q only in the case where we want an estimation of properties for cloth reinforcement
        it has to use the micromechanical model in calc_local_properties"""
        if self.Qvalid:
            return self.Q
        properties = self.calc_local_properties()
        properties['poisson21'] = self.calc_poisson21()
        Q = np.zeros((3, 3))
        Q[0, 0] = properties['E11'] / (1 - properties['poisson12'] * properties['poisson21'])
        Q[0, 1] = Q[1, 0] = (properties['poisson21'] * properties['E11'] /
                             (1 - properties['poisson12'] * properties['poisson21']))
        Q[1, 1] = properties['E22'] / (1 - properties['poisson12'] * properties['poisson21'])
        a = 1.0
        if self.arch == 'satin' and len(self.angles) == 2:
            a = 1.2
        elif self.arch == 'twill' and len(self.angles) == 2:
            a = 1.5
        elif self.arch == 'plain' and len(self.angles) == 2:
            a = 2
        elif len(self.angles) != 2 and self.arch != 'uni':
            print('There seem to be an error in the specification of reinforcement architecture')
            print('For one angle only check that architecture is unidirectional')
            print('For textile reinforcement check that two and only two angles are specified')
            raise
        Q[2, 2] = a * properties['G12']
        return Q

    def calc_Q(self):
        """Since we need calc_micro_Q_matrix to be in this class and Q calculation requires access to lamina
        properties we also put here this function"""
        Q = np.zeros((3, 3))
        properties = {}
        for property in self.xml:
            try:
                properties[property.tag] = float(property.text)
            except:
                pass

        Q[0, 0] = properties['E11'] / (1 - properties['poisson12'] * properties['poisson21'])
        Q[0, 1] = Q[1, 0] = (properties['poisson21'] * properties['E11'] /
                             (1 - properties['poisson12'] * properties['poisson21']))
        Q[1, 1] = properties['E22'] / (1 - properties['poisson12'] * properties['poisson21'])
        Q[2, 2] = properties['G12']
        self.Q = Q  # Update Q and keep it to avoid recalculation in a multiply laminate
        self.Qvalid = True
        return Q

    def T_matrix(self, theta):
        """used to calculate the properties of a lamina once it has been rotated. It is here instead of in ply
        to be used to calculate properties of multiaxial cloth reinforcement"""
        th = np.deg2rad(theta)
        T = np.zeros((3, 3))
        T[0, 0] = np.cos(th)**2
        T[0, 1] = np.sin(th)**2
        T[0, 2] = 2 * np.sin(th) * np.cos(th)
        T[1, 0] = np.sin(th)**2
        T[1, 1] = np.cos(th)**2
        T[1, 2] = - 2 * np.sin(th) * np.cos(th)
        T[2, 0] = - np.sin(th) * np.cos(th)
        T[2, 1] = np.sin(th) * np.cos(th)
        T[2, 2] = np.cos(th)**2 - np.sin(th)**2
        return T

    def calc_multiaxial_properties(self):
        A = np.zeros((3, 3))
        Q = self.calc_micro_Q_matrix()
        #TODO formular diferente para eliminar la dependencia de t
        t = self.t
        for angle in self.angles:
            T = self.T_matrix(-angle)
            Qbar = T.dot(Q).dot(T.T)
            A += t * Qbar
        Abar = np.linalg.inv(A)
        # TODO formular diferente para eliminar la dependencia de t
        total_t = t * len(self.angles)
        E11 = 1 / (total_t * Abar[0, 0])
        E22= 1 / (total_t * Abar[1, 1])
        G12 = 1 / (total_t * Abar[2, 2])
        poisson12 = - Abar[0, 1] / Abar[0, 0]
        poisson21 = poisson12 * E22 / E11
        return {'E11': E11, 'E22': E22, 'G12': G12, 'poisson12': poisson12, 'poisson21': poisson21}

if __name__ == '__main__':
    l = Lamina()
    l.set_micro_prop(arch='satin', angles=[0,90])
    print(l.calc_micro_Q_matrix())
    l.calc_Q()