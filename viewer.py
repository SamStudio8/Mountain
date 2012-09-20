#!/usr/bin/env python
"""A Python Curses genome viewer.
Displays each sample in a directory of FASTA samples or a FASTA library on
each line and allows traversal of the sequences."""

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
#
#This program is distributed in the hope that it will be useful, but WITHOUT
#ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
#details.
#
#You should have received a copy of the GNU Lesser General Public License along
#with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import os
import re
import curses

def execute(stdscr, args):
    NAME_SIZE = 20
    CONVERT_LINES = re.compile("\r?\n|\r")
    curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_WHITE)
    curses.init_pair(2, curses.COLOR_WHITE, curses.COLOR_RED)

    ref = open(args.reference, 'r')
    ref_readname = ref.readline().strip()
    ref_seq = CONVERT_LINES.sub("", ref.read())
    largest_string = len(ref_seq)

    if args.folder:
        fasta_list = sorted(os.listdir(args.folder))
    
    snps = []
    if args.snps:
        handler = open(args.snps, 'r')
        for line in handler:
            try:
                snps.append(int(line.strip()))
            except:
                pass

    sequences = []
    curr_iter = 0
    for fasta in fasta_list:
        curr_iter += 1
        f = {}
        handler = open(args.folder+fasta, 'r')
        readname = handler.readline().strip()
        stdscr.clear() 
        stdscr.addstr(0,0,"[LOAD]\tStand Back. I am a loadings...")
        stdscr.addstr(1,0,"[READ]["+str(curr_iter)+" of "+str(len(fasta_list))+"]\t"+readname)
        stdscr.refresh()
        seq = CONVERT_LINES.sub("", handler.read())
        f['readname'] = readname
        f['sequence'] = seq
        sequences.append(f)

        if len(seq) > largest_string:
            largest_string = len(seq)
    
    offset = 0
    vertical_offset = 0
    while True:
        stdscr.clear() 
        stdscr.addstr(0, 0, "POS")
        stdscr.addstr(1, 0, ref_readname[0:19])
        stdscr.addstr(1, NAME_SIZE, ref_seq[offset:(offset+stdscr.getmaxyx()[1]-NAME_SIZE)], curses.color_pair(1))
        current_line = 2

        for sample in sequences[vertical_offset : (vertical_offset + stdscr.getmaxyx()[0])]:
            try:
                stdscr.addstr(current_line, 0, sample['readname'][0:19])
                for index, base in enumerate(sample['sequence'][offset:(offset+stdscr.getmaxyx()[1]-NAME_SIZE)]):
                    if (offset+index+1) in snps:
                        #TODO Addchar?
                        stdscr.addstr(current_line, NAME_SIZE+index, base, curses.color_pair(2))
                    else:
                        if base != ref_seq[index+offset]:
                            stdscr.addstr(current_line, NAME_SIZE+index, base, curses.color_pair(2))
                        else:
                            stdscr.addstr(current_line, NAME_SIZE+index, base)
                current_line += 1
            except:
                pass

        event = stdscr.getch() 
        if event == ord("q"): break 
        elif event == curses.KEY_LEFT:
            if offset > 0:
                offset += -1

        elif event == curses.KEY_RIGHT: 
            if offset < (largest_string - stdscr.getmaxyx()[1] + NAME_SIZE):
                offset += 1

        elif event == curses.KEY_UP:
            if vertical_offset > 0:
                vertical_offset += -1

        elif event == curses.KEY_DOWN:
            if vertical_offset < len(sequences):
                vertical_offset += 1
                
        elif event == curses.KEY_HOME:
            offset = 0

        elif event == curses.KEY_END:
            offset = (largest_string - stdscr.getmaxyx()[1] + NAME_SIZE)

        # Next Block
        elif event == curses.KEY_PPAGE:
            if ((largest_string - (offset + stdscr.getmaxyx()[1])) > 
                    (stdscr.getmaxyx()[1] - NAME_SIZE)):
                offset += (stdscr.getmaxyx()[1] - NAME_SIZE)
            else:
                offset = (largest_string - stdscr.getmaxyx()[1] + NAME_SIZE)

        # Prev Block
        elif event == curses.KEY_NPAGE:
            if offset > (stdscr.getmaxyx()[1]):
                offset += -(stdscr.getmaxyx()[1] - NAME_SIZE)
            else:
                offset = 0

    curses.nocbreak(); stdscr.keypad(0); curses.echo()
    curses.endwin()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(                                              
        description="Output series of subsequences for manual comparison.")
    parser.add_argument('folder', metavar='folder', type=str,                      
        help="Folder of FASTA sequences to compare")                               
    parser.add_argument('reference', metavar='reference', type=str,                
        help="FASTA file to use as reference")                                     
    parser.add_argument('--snps', metavar="", type=str)
    #TODO Support library files (abstract FASTA reading!)
    #TODO Optional "SNP Locations" file
    args = parser.parse_args()                                                
    curses.wrapper(execute, args)
