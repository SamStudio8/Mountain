#!/usr/bin/env python
"""A set of classes to enable the handling of FASTA files and libraries"""

__author__ = "Sam Nicholls <sn8@sanger.ac.uk>"
__copyright__ = "Copyright (c) 2012 Genome Research Ltd."
__version__ = 0.15
__license__ = "GNU Lesser General Public License V3"
__maintainer__ = "Sam Nicholls <sam@samnicholls.net>"

#This file is part of Mountain.
#
#Mountain is free software: you can redistribute it and/or modify it under 
#the terms of the GNU Lesser General Public License as published by the Free 
#Software Foundation; either version 3 of the License, or (at your option) any  
#later version.

#This program is distributed in the hope that it will be useful, but WITHOUT 
#ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
#FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more 
#details.

#You should have received a copy of the GNU Lesser General Public License along 
#with this program.  If not, see <http://www.gnu.org/licenses/>.

import re

CONVERT_LINES = re.compile("\r?\n|\r")

class FASTA_Library(object):
    """Wraps a file handler and provides iterable access over each FASTA"""

    def __init__(self, filepath):
        """Constructs the read only file handler and parses the header"""
        self.handler = open(filepath, 'r')
        self.process_file()

    def process_file(self):
        """Ensure the file contains FASTA sequences and remove any headers"""
        self.handler.seek(0)
        line = self.handler.readline().strip()

        # Remove any possible headers
        while(line[0] != '>'):
            line = self.handler.readline()

            if not line:
                print ("Reached EOF without finding a readname. "
                    "No FASTA sequences found.")
                exit(0)

        self.current_readname = line
        self.has_next = 1

    def close(self):
        """Close the file handler"""
        self.handler.close()

    def __iter__(self):
        self.process_file() #Reset file pointer and re-process file header
        return self

    def next(self):
        if self.has_next:
            name = self.current_readname.strip()
            self.has_next = 0
        else:
            raise StopIteration

        seq = ""
        line = self.handler.readline()
        while line and self.has_next == 0:
            if line[0] != '>':
                seq += CONVERT_LINES.sub("", line)
                line = self.handler.readline()
            else:
                self.current_readname = line
                self.has_next = 1

        return FASTA(name, seq) 

    def get_next(self):
        """Return the next FASTA sequence in the library, if there is one.
        Typically used in place of an iterator where you know there is only
        one sample inside the file"""
        if self.has_next:
            return self.next()
        else:
            return None

class FASTA(object):
    """Simple encapsulation of a FASTA sequence"""

    def __init__(self, name, seq):
        """Construct the FASTA with a readname and sequence"""
        self.readname = name
        self.seq = seq

    def __iter__(self):
        """Provide an iterator over the bases in this FASTA sequence"""
        return iter(self.seq)

    def __len__(self):
        """Return the length of this FASTA seqeunce"""
        return len(self.seq)
