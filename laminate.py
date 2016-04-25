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
from ply import Ply
from lamina import Lamina
import xml.etree.ElementTree as ET


class Laminate:
    """Class to hold laminate definition. To help with parsing it to/from a file we use an xml tree to hold
    its data. This also comes from using elamX as a GUI. Since data is already in the xml tree we use it.
    """
    def __init__(self, name='Laminate name', uuid='laminate identifier'):
        self.xml = ET.Element('laminate')
        self.xml.set('name', name)
        self.xml.set('uuid', uuid)
        self.plies = []  # List of ply objects that form the laminate
        self.S_bars = []  # List for convenience of calculating engineering properties
        self.ts = []  # List for convenience of calculating engineering properties
        self.A = np.zeros((3, 3))
        self.B = np.zeros((3, 3))
        self.D = np.zeros((3, 3))
        self.symmetric = False
        self.middle_layer = False
        self.Ex = 0.0
        self.Ey = 0.0
        self.Gxy = 0.0
        self.poissonxy = 0.0
        self.zx = 0.0
        self.zy = 0.0
        self.EIx = 0.0
        self.EIy = 0.0
        self.GA = 0.0
        self.total_t = 0.0

    def calc_laminate(self):
        self.A = 0.0
        self.B = 0.0
        self.D = 0.0
        self.ts = []
        self.S_bars = []
        total_lam = total_lam = self.plies
        if self.symmetric:
            if self.middle_layer:
                total_lam = self.plies + list(reversed(self.plies))[1:]
            else:
                total_lam = self.plies + list(reversed(self.plies))
        for ply in total_lam:
            self.ts.append(float(ply.xml.find('thickness').text))
        h = -sum(self.ts) / 2.0  # Distante from bottom of surface to middle surface following Barbero's
        for ply in total_lam:
            t = float(ply.xml.find('thickness').text)
            self.S_bars.append(ply.calc_S_bar())
            Q_bar = ply.calc_Q_bar()
            self.A += ((h + t) - h) * Q_bar
            self.B += 1 / 2.0 * ((h + t)**2 - h**2) * Q_bar
            self.D += 1 / 3.0 * ((h + t)**3 - h**3) * Q_bar
            h += t
        self.total_t = 2 * h
        self.calc_enginprops()
        return

    def calc_enginprops(self):
        """ Given matrices A, B, and D calculates equivalent engineering properties.
        It's separated from matix calculations for convenience.
        Calling it by itself is not needed """
        # Let's assemble the ABD matrix even if it is not required
        ABD = np.bmat([[self.A, self.B], [self.B, self.D]])
        ABD_inv = np.linalg.inv(ABD)
        # We would use the whole matrix. This gives results similar to elamX and considers poisson effects
        A_inv = ABD_inv[0:3, 0:3]
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
            Ex = 1 / S_bar[0, 0]
            Ey = 1 / S_bar[1, 1]
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
            Ex = 1 / S_bar[0, 0]
            Ey = 1 / S_bar[1, 1]
            Gxy = 1 / S_bar[2, 2]
            z += t / 2.0
            self.EIx += Ex * (t**3 / 12 + t * (z - self.zx)**2)
            self.EIy += Ey * (t**3 / 12 + t * (z - self.zy)**2)
            self.GA += Gxy * t
            z += t / 2.0
        return self.Ex, self.Ey, self.Gxy, self.poissonxy


def parse_elamx_file(filename):
    """Utility function to parse a elamx file into a list of objects of type laminate"""
    laminate_list = []
    lamina_dict = {}
    tree = ET.parse(filename)
    root = tree.getroot()
    laminates = root.find('laminates')  # laminates in the file
    materials = root.find('materials')  # laminae in the file
    for material in materials:
        lamina_dict[material.attrib['uuid']] = parse_elamx_material(material)
    for laminate in laminates:
        laminate_list.append(parse_elamx_laminate(laminate, lamina_dict))
    return laminate_list


def parse_elamx_laminate(laminate_elamx, lamina_dict):
    name = laminate_elamx.attrib['name']
    uuid = laminate_elamx.attrib['uuid']
    temp_laminate = Laminate(name=name, uuid=uuid)
    if laminate_elamx.attrib['symmetric'] == 'true':
        temp_laminate.symmetric = True
    if laminate_elamx.attrib['with_middle_layer'] == 'true':
        temp_laminate.middle_layer = True
    for layer in laminate_elamx:
        t = float(layer.find('thickness').text)
        angle = float(layer.find('angle').text)
        nam = layer.attrib['name']
        uuid = layer.attrib['uuid']
        mat = layer.find('material').text
        lamina = lamina_dict[mat]
        temp_ply = Ply(lamina, name=nam, uuid=uuid, t=t, angle=angle)
        temp_laminate.plies.append(temp_ply)
        temp_laminate.xml.append(lamina.xml)
    return temp_laminate


def parse_elamx_material(material):
    name = material.attrib['name']
    uuid = material.attrib['uuid']
    temp_lamina = Lamina(name=name, uuid=uuid)
    properties = ['E11', 'E22', 'G12', 'poisson12', 'rho', 'sigma1_t', 'sigma1_c',
                  'sigma2_t', 'sigma2_c', 'tau12']
    elamx_prop = ['Epar', 'Enor', 'G', 'nue12', 'rho', 'RParTen', 'RParCom',
                  'RNorTen', 'RNorCom', 'RShear']
    for prop, eprop in zip(properties, elamx_prop):
        temp_lamina.xml.find(prop).text = material.find(eprop).text
    #  Set poisson 21 as calculated value
    temp_lamina.xml.find('poisson21').text = str(temp_lamina.calc_poisson21())
    return temp_lamina


def read_laminate_by_name(elamx_file, lam_name):
    """Function to read a laminate from the file given its name.
    No effort is done to check if more than one laminate have the same name.
    In the file only uuid are unique"""
    lams = parse_elamx_file(elamx_file)
    for laminate in lams:
        if laminate.xml.attrib['name'] == lam_name:
            return laminate
    return

if __name__ == '__main__':
    print('Testing laminate reader')
    all_lams = parse_elamx_file('example.elamx')
    print('Names of detected laminates:')
    for l in all_lams:
        print(l.xml.attrib['name'])
    print('uuid of detected laminates:')
    for l in all_lams:
        print(l.xml.attrib['uuid'])
    print('Finding a particular laminate')
    if read_laminate_by_name('example.elamx', 'Laminate2'):
        print('Laminate found')
    test_lam = read_laminate_by_name('example.elamx', 'Barbero')
    test_lam.calc_laminate()
    print(test_lam.A)
    print(test_lam.B)
    print(test_lam.D)
    print('Ex: ' + str(test_lam.Ex))
    print('Ey: ' + str(test_lam.Ey))
    print('Gxy: ' + str(test_lam.Gxy))
    print('poissonxy: ' + str(test_lam.poissonxy))
    print('EIx: ' + str(test_lam.EIx))
