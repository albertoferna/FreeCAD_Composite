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

class Ply:
    """This class hold a lamina of certain thickness oriented in a certain direction. It takes its data from
    a lamina and calculates the needed matrices."""

    def __init__(self, lamina, name='Layer', uuid='', t=0.1, angle=0.0):
        self.xml = ET.Element('ply')
        self.set('name', name)
        self.set('uuid', uuid)
        self.lamina = lamina

        self.xml.append('thickness')
        self.xml.find('thickness').text = str(t)
        self.xml.append('angle')
        self.xml.find('angle').text = str(angle)
        self.xml.append('lamina_id')
        self.xml.find('lamina_id').text = lamina.xml.attrib['uuid']



