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

from lamina import Lamina
import xml.etree.ElementTree as ET
import numpy as np


class Ply:
    """This class hold a lamina of certain thickness oriented in a certain direction. It takes its data from
    a lamina and calculates the needed matrices."""

    def __init__(self, lamina, name='Layer', uuid='', t=0.1, angle=0.0):
        self.xml = ET.Element('ply')
        self.xml.set('name', name)
        self.xml.set('uuid', uuid)
        self.lamina = lamina
        self.xml.append(ET.Element('thickness'))
        self.xml.find('thickness').text = str(t)
        self.xml.append(ET.Element('angle'))
        self.xml.find('angle').text = str(angle)
        self.xml.append(ET.Element('lamina_id'))
        self.xml.find('lamina_id').text = lamina.xml.attrib['uuid']
        self.Q_bar = np.zeros((3, 3))

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

    def calc_Q_bar(self):
        Q = self.lamina.calc_Q()
        theta = float(self.xml.find('angle').text)
        T_inv = self.calc_T_matrix(-theta)  # as a rotation the inverse is the inverse rotation
        Q_bar = T_inv.dot(Q).dot(T_inv.T)  # Instead of using Reuter's matrix we can use (T**-1).transpose (Barbero)
        return Q_bar

    def calc_S_bar(self):
        S_bar = np.linalg.inv(self.calc_Q_bar())
        return S_bar


if __name__ == '__main__':
    l = Lamina(uuid='lamina identifier')
    l.set_micro_prop()
    p = Ply(l, name='first ply', uuid='ply identifier', t=0.5, angle=45)
    ET.dump(p.xml)
    print(p.calc_S_bar())
