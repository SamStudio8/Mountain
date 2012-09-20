#!/usr/bin/env python
"""A set of classes to enable the handling of aligned block data from psl and 
pslx files."""

__author__ = "Sam Nicholls <sn8@sanger.ac.uk>"
__copyright__ = "Copyright (c) 2012 Genome Research Ltd."
__version__ = 0.50
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

class PSL(object):
    """Wraps a file handler and provides iterable access over each 'BlockSet' 
    (list of successfully aligned Blocks) in the file"""

    def __init__(self, filepath):
        """Constructs the read only file handler and parses the header"""
        self.handler = open(filepath, 'r')
        self.handler.seek(0)
        line = self.handler.readline().strip()
        if line == "psLayout version 3":
            while(line[0] != '-'):
                line = self.handler.readline()
        else:
            print "Unknown PSL format! Expecting 'psLayout version 3'"
            exit(0)

    def close(self):
        """Close the file handler"""
        self.handler.close()

    def __iter__(self):
        self.handler.seek(0) # Reset the file pointer and re-process header
        return self

    def next(self):
        fields = self.handler.readline().strip().split("\t")
        if len(fields) >= 21:
            if len(fields) > 21:
                block_q_seqs_i = fields[21].upper().split(',')[:-1]
                block_t_seqs_i = fields[22].upper().split(',')[:-1]
            else:
                block_q_seqs_i = [''] * int(fields[17])
                block_t_seqs_i = [''] * int(fields[17])

            return BlockSet(
                match = int(fields[0]),
                mismatch = int(fields[1]),
                repmatch = int(fields[2]),
                n = int(fields[3]),
                q_gap_count = int(fields[4]),
                q_gap_bases = int(fields[5]),
                t_gap_count = int(fields[6]),
                t_gap_bases = int(fields[7]),
                strand = fields[8],
                q_name = fields[9],
                q_size = int(fields[10]),
                q_start = int(fields[11]),
                q_end = int(fields[12]),
                t_name = fields[13],
                t_size = int(fields[14]),
                t_start = int(fields[15]),
                t_end = int(fields[16]),
                block_count = int(fields[17]),
                block_sizes = map(int, fields[18].split(',')[:-1]),
                block_q_starts = map(int, fields[19].split(',')[:-1]),
                block_t_starts = map(int, fields[20].split(',')[:-1]),
                block_q_seqs = block_q_seqs_i,
                block_t_seqs = block_t_seqs_i)
        else:
            raise StopIteration
        

class BlockSet(object):
    """Holds the set of query subsequences and their mappings to the target
    reference sequence for a particular sample. Provides iterable access over 
    each Block in this BlockSet"""

    def __init__(self, **kwds):
        """Construct the BlockSet"""
        #TODO Should provide a way to check all required variables were
        # assigned to ensure future readers create BlockSets of the same form.
        self.__dict__.update(kwds)

    def __iter__(self):
        self.__iter_pointer = -1
        return self

    def next(self):
        self.__iter_pointer += 1
        if self.__iter_pointer < self.block_count:
            return self.get_block(self.__iter_pointer)
        else:
            raise StopIteration

    def get_block(self, i):
        """Get the ith Block in this BlockSet, if it exists"""
        if i < self.block_count:
            return Block(
                i,
                self.block_sizes[i],
                self.block_q_starts[i],
                self.block_t_starts[i],
                self.block_q_seqs[i],
                self.block_t_seqs[i])
        else:
            return None

    def get_all_blocks(self):
        """Get all Blocks in this BlockSet"""
        blocks = []
        for i in range(0, self.block_count):
            blocks.append(self.get_block(i))
        return blocks

class Block(object):
    """An individual alignment; representing a subsequence of the query
    sequence and how it maps to the target reference sequence"""

    def __init__(self, i, size, q_st, t_st, q_seq, t_seq):
        """Construct the Block with its position in the blockset, size, 
        1-index query sequence start, 1-index target sequence start and if the
        read came from a pslx, the aligned query and target subsequences"""
        self.no = i
        self.size = size
        self.q_start = q_st
        self.t_start = t_st
        self.q_seq = q_seq
        self.t_seq = t_seq

    def is_first(self):
        """Check if this Block is the first in the BlockSet"""
        return self.no == 0

    def is_last(self, blockset_len):
        """Check if this Block is the last in the Blockset given the length 
        of the BlockSet"""
        return self.no == (blockset_len - 1)

