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

import numpy as np



class Laminate:
    """Class to hold laminate definition. It contains a xml tree that can be read from elamx
    """

    materials = {}  # keeps a dict of used materials
    S_bars = []  # keeps rotated S matrices to save computing time
    ts = []  # thickness also kept here for loop convenience
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    symmetric = False
    sandwich = ''
    Ecore = 0.0
    Gcore = 0.0
    Ex = 0.0
    Ey = 0.0
    Gxy = 0.0
    poissonxy = 0.0
    total_t = 0.0
    t_core = 0.0
    zx = 0.0
    zy = 0.0
    EIx = 0.0
    EIy = 0.0
    GA = 0.0

    def read_materials(self, elamx):
        """Reads materials present in an xml structure that has been read from an elamx file"""
        return
    def insert_lamina(self, lamina, angle, pos):
        """ Takes a lamina and adds its properties to the corresponding list
        :param Lamina, angle, position
        :return: None
        """
        self.layers.insert(pos, (lamina, angle))

    def del_lamina(self, pos):
        """ Takes a lamina and adds its properties to the corresponding list
        :param Lamina, angle, position
        :return: None
        """
        del self.layers[pos]

    def read_laminate(self, elamx_file, lam_uuid):
        """Function to read a laminate from the file given its uuid. In the file only uuid are unique"""
        self.layers = []
        self.materials = {}
        self.symmetric = False
        f = open(elamx_file, 'r')
        soup_xml = BeautifulSoup(f, 'lxml')
        laminate = soup_xml.find('laminate', attrs={'uuid': lam_uuid})
        self.symmetric = laminate.attrs['symmetric'] == 'true'
        for layer in laminate.findAll('layer'):
            self.layers.append(layer)
            mat_id = layer.material.text
            self.materials[mat_id] = (soup_xml.find(uuid=mat_id))
            # check for core material
            coreCheck = layer.attrs['name'].lower()
            if 'honeycomb' in coreCheck:
                self.sandwich = 'honeycomb'
                self.Ecore = float(self.materials[mat_id].epar.text)
                self.t_core = float(layer.thickness.text)
                self.Gcore = float(self.materials[mat_id].g.text)
            elif 'foam' in coreCheck:
                self.sandwich = 'foam'
                self.Ecore = float(self.materials[mat_id].epar.text)
                self.t_core = float(layer.thickness.text)
                self.Gcore = float(self.materials[mat_id].g.text)
        if self.symmetric:
            for layer in reversed(laminate.findAll('layer')):
                self.layers.append(layer)
                mat_id = layer.material.text
                self.materials[mat_id] = (soup_xml.find(uuid=mat_id))
        return

    def read_laminate_by_name(self, elamx_file, lam_name):
        """Function to read a laminate from the file given its name.
        No effort is done to check if more than one laminate have the same name.
        In the file only uuid are unique"""
        f = open(elamx_file, 'r')
        soup_xml = BeautifulSoup(f, 'lxml')
        laminate = soup_xml.find('laminate', attrs={'name': lam_name})
        uuid = laminate.attrs['uuid']
        self.read_laminate(elamx_file, uuid)
        return

    def calc_laminate_elamx(self):
        self.A = 0.0
        self.B = 0.0
        self.D = 0.0
        self.S_bars = []
        Qs = []
        angles = []
        self.ts = []
        for layer in self.layers:
            mat_id = layer.material.text
            E11 = float(self.materials[mat_id].epar.text)
            E22 = float(self.materials[mat_id].enor.text)
            G12 = float(self.materials[mat_id].g.text)
            poisson12 = float(self.materials[mat_id].nue12.text)
            Qs.append(self.calc_Q(E11, E22, G12, poisson12))
            angles.append(float(layer.angle.text))
            self.ts.append(float(layer.thickness.text))
        h = - sum(self.ts) / 2.0
        for Q, angle, t in zip(Qs, angles, self.ts):
            T_inv = self.calc_T_matrix(-angle)  # as a rotation the inverse is the inverse rotation
            Q_bar = T_inv.dot(Q).dot(T_inv.T)  # Instead of using Reuter's matrix we can use (T**-1).transpose (Barbero)
            self.S_bars.append(np.linalg.inv(Q_bar))  # Save calculation by doing it here and storing them
            self.A += ((h + t) - h) * Q_bar
            self.B += 1 / 2 * ((h + t)**2 - h**2) * Q_bar
            self.D += 1 / 3 * ((h + t)**3 - h**3) * Q_bar
            h += t
        self.total_t = 2 * h
        self.calc_enginprops()
        return

    def calc_laminate(self):
        self.A = 0.0
        self.B = 0.0
        self.D = 0.0
        self.S_bars = []
        Qs = []
        angles = []
        self.ts = []
        for layer in self.layers:
            mat_id = layer[0].material_id
            E11 = layer[0].E11
            E22 = layer[0].E22
            G12 = layer[0].G12
            poisson12 = layer[0].poisson12
            Qs.append(self.calc_Q(E11, E22, G12, poisson12))
            angles.append(layer[1])
            self.ts.append(layer[0].t)
        h = - sum(self.ts) / 2.0
        for Q, angle, t in zip(Qs, angles, self.ts):
            T_inv = self.calc_T_matrix(-angle)  # as a rotation the inverse is the inverse rotation
            Q_bar = T_inv.dot(Q).dot(T_inv.T)  # Instead of using Reuter's matrix we can use (T**-1).transpose (Barbero)
            self.S_bars.append(np.linalg.inv(Q_bar))  # Save calculation by doing it here and storing them
            self.A += ((h + t) - h) * Q_bar
            self.B += 1 / 2 * ((h + t)**2 - h**2) * Q_bar
            self.D += 1 / 3 * ((h + t)**3 - h**3) * Q_bar
            h += t
        self.total_t = 2 * h
        self.calc_enginprops()
        return

    def calc_enginprops(self):
        """ Given matrices A, B, and D calculates equivalent engineering properties.
        It's separated from matix calculations for convenience.
        Calling it by itself is not needed """
        # Let's assemble the ABD matrix even if it is not required in the guide
        ABD = np.bmat([[self.A, self.B], [self.B, self.D]])
        ABD_inv = np.linalg.inv(ABD)
        # As required by the guide we would not use matrix B.
        # This is equivalent to have a symmetric laminate with the same total thickness
        A_inv = np.linalg.inv(self.A)
        self.Ex = 1 / (self.total_t * A_inv[0, 0])  # It is 2 * t because we need total thickness
        self.Ey = 1 / (self.total_t * A_inv[1, 1])
        self.Gxy = 1 / (self.total_t * A_inv[2, 2])
        self.poissonxy = - A_inv[0,1] / A_inv[0, 0]
        # Flexural stiffness properties
        self.zx = 0.0
        self.zy = 0.0
        zx_dem = 0.0
        zy_dem = 0.0
        self.EIx = 0.0
        self.EIy = 0.0
        z = 0.0
        # Calculate neutral axis in direction x and y
        for S_bar, t in zip(self.S_bars, self.ts):
            Ex = 1 / S_bar[0,0]
            Ey = 1 / S_bar[1,1]
            z += t / 2.0
            self.zx += Ex * t * z
            zx_dem += Ex * t
            self.zy += Ey * t * z
            zy_dem += Ey * t
            z += t / 2.0
        self.zx = self.zx / zx_dem
        self.zy = self.zy / zy_dem
        # Calculate EI in direction x and y
        z = 0.0
        for S_bar, t in zip(self.S_bars, self.ts):
            Ex = 1 / S_bar[0,0]
            Ey = 1 / S_bar[1,1]
            Gxy = 1 / S_bar[2,2]
            z += t / 2.0
            self.EIx += Ex * (t**3 / 12 + t * (z - self.zx)**2)
            self.EIy += Ey * (t**3 / 12 + t * (z - self.zy)**2)
            self.GA += Gxy * t
            z += t / 2.0
        return self.Ex, self.Ey, self.Gxy, self.poissonxy

    def calc_Q(self, E11, E22, G12, poisson12):
        Q = np.zeros((3, 3))
        poisson21 = poisson12 * E22 / E11
        Q[0, 0] = E11 / (1 - poisson12 * poisson21)
        Q[0, 1] = Q[1, 0] = poisson21 * E11 / (1 - poisson12 * poisson21)
        Q[1, 1] = E22 / (1 - poisson12 * poisson21)
        Q[2, 2] = G12
        return Q

    def calc_T_matrix(self, theta):
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
