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

from lamina import Lamina
import xml.etree.ElementTree as ET

class Lamina3D(Lamina):
    """This class serves a function very similar to lamina except it is used to calculate the full compliance and
    stiffness matrix. It follows the procedure in "Three-Dimensional Effective Property and Strength Prediction of
    Thick Laminated Composite Media" by Bogetti, Hoppel, and Drysdale from the army research laboratory.
    """
    def __init__(self):
        Lamina.__init__(self, name="nombre de prueba")
        self.xml.append(ET.Element("poisson23"))

    def calc_C(self):
        """Calculates stiffness matrix for a transversely isotropic lamina. It takes all the properties from
        a general lamina description. The remaining value poisson23 is assumed to be 1.65 time poisson12"""
        pass
if __name__ == '__main__':
    l = Lamina3D()
    l.set_micro_prop()
    print(l.xml.get(('name')))
    print(l.xml.find('E11').text)

