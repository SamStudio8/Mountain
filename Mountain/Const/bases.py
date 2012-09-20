#!/usr/bin/env python
"""Standard Bases and IUPAC Lookup Dictionary"""

__author__ = "Sam Nicholls <sn8@sanger.ac.uk>"
__copyright__ = "Copyright (c) 2012 Genome Research Ltd."
__version__ = 1.0
__license__ = "GNU Lesser General Public License V3"
__maintainer__ = "Sam Nicholls <sam@samnicholls.net>"

#This file is part of Mountain.
#
#Mountain is free software: you can redistribute it and/or modify it under
#the terms of the GNU Lesser General Public License as published by the Free
#Software Foundation; either version 3 of the License, or (at your option) any
#later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT
#ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
#details.
#
#You should have received a copy of the GNU Lesser General Public License along
#with this program.  If not, see <http://www.gnu.org/licenses/>.

# Standard Bases
STANDARD_BASES = ['G', 'C', 'T', 'A']

# Lookup dict for reconciling IUPAC codes to standard bases.
# Constructed with data from:
# http://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation
IUPAC = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'T', 'G']}
