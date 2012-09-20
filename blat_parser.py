#!/usr/bin/env python
"""Parses alignment output from blat. For each line the set of aligned block 
data is applied to the query sequence to construct a new reference aligned 
sequence. Query sequence data can be read from a directory of extracted 
FASTA files, a FASTA library or the subsequence output from a pslx."""

__author__ = "Sam Nicholls <sn8@sanger.ac.uk>"
__copyright__ = "Copyright (c) 2012 Genome Research Ltd."
__version__ = 0.975
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

import textwrap
import re
import os
import argparse

from Mountain.IO import blatreaders, fastareaders

def execute(args):
    convertLines = re.compile("\r?\n|\r")
    """Execute the script with argpased arguments"""

    # If not using subsequences from a pslx file, build a FASTA readname to
    # FASTA filename mapping for all files in the provided list to allow 
    # script to reconcile the sample names to their actual FASTA later
    files = {}
    if not args.pslx:
        print "[READ]\tImporting Extracted FASTA File List"
        filelist = open(args.FASTA_list, 'r')
        for f in filelist:
            current = open(args.FASTA_path+f.strip(), 'r')
            files[current.readline().strip().split(' ')[0]] = f.strip()
            current.close()

    readname_list = []
    blocksets = {}

    # Load the file reader, currently only supports psl/pslx
    if args.type == "psl" or args.type == "pslx" or args.type == "psl/pslx":
        blat_output = blatreaders.PSL(args.blat_out)
    else:
        print "Filetype not supported. Use --type[psl|pslx]."
        exit(0)

    print "[READ]\tParsing BLAT Output"
    for blockset in blat_output:
        if blockset.strand == '-':
            #TODO Complement handling
            continue

        if args.multi:
            blockset_name = blockset.q_name+"::"+blockset.t_name
        else:
            blockset_name = blockset.q_name

        if blockset_name not in readname_list:
            #TODO Multistrand handling
            readname_list.append(blockset_name)
            blocksets[blockset_name] = blockset

    print "[OKAY]\t"+str(len(blocksets))+" Blocksets Read"
    print "[READ]\tProcessing FASTA files"

    loopPos = 0
    for blockset_name in blocksets:
        bset = blocksets[blockset_name]
        if not args.pslx:
            if '>'+blockset_name not in files:
                continue

            loopPos += 1
            fasta_lib = fastareaders.FASTA_Library(
                args.FASTA_path+files['>'+blockset_name])
            fasta = fasta_lib.get_next()
            name = fasta.readname
        else:
            name = '>'+blockset_name
            loopPos += 1

        termination = 0
        
        if not os.path.isdir(args.aligned):
            os.makedirs(args.aligned)
            
        if not args.merged:
            if not args.pslx:
                output = open(args.aligned+files['>'+blockset_name], 'w+')
            else:
                output = open(args.aligned+blockset_name+".fa", "w+")
            output.write(name+"\n")
        else:
            try:
                output.write(">"+blockset_name+"\n")
            except NameError: #TODO Terrible
                #TODO Allow entering filename!
                output = open(args.aligned+"output", "w+")
                output.write(">"+blockset_name+"\n")

        # Iterate over each BLAT block for this particular sequence
        pieces = []
        for block in bset:
            # Extract the current block from the current query FASTA
            if not args.pslx:
                piece = list(fasta.seq[block.q_start : block.q_start+block.size])
            else:
                piece = list(block.q_seq)

            # If this is the first block and the query sequence wasn't mapped
            # to the start of the target reference, pad the leading gap with
            # the padding character (_).
            #
            # Q |==========>
            #      \\\\\\\\
            # T |XX========>
            #
            if block.is_first() and block.t_start != 0:
                piece[0] = '*'
                piece = list(block.t_start * '_')+piece


            # If the previous block completed with a termination, mark the
            # start of this block with the termination character (*) also.
            if(termination==1):
               piece[0] = '*'
               termination = 0

            # If this is not the last block, check for insertions and
            # deletions between this block and the next.
            # Note that an insertion and deletion can occur at the same time!
            if not block.is_last(bset.block_count):

                # Insertion: Base(s) in query not found in target (ref).
                #
                # Q >====O======>
                #    ||||//////
                # T >==========>
                #
                # If the query sequence for the next block doesn't begin 
                # immediately following the position where this query block
                # ends, the bases between the two positions were not mapped to 
                # the target, and are thus insertions.
                # These are not appended to the new aligned sequence and the
                # end of this block is marked with the termination character.
                if( (block.q_start + block.size) 
                        < bset.get_block(block.no+1).q_start ):
                    piece[-1] = '*'
                    termination = 1

                # Deletion: Base(s) in target (ref) not found in query.
                #
                # Q >=====>
                #    ||\\\
                # T >==X===>
                #
                # If the target sequence for the next block doesn't begin
                # immediately following the position where this target 
                # sequence block ends, the bases between the two points were
                # not mapped from the query and are thus deletions.
                # To maintain alignment, the "missing" bases are appended to 
                # the end of this block with the padding character (_).
                # The last character will be the termination character (*).
                if ( (block.t_start + block.size)
                        < bset.get_block(block.no+1).t_start ):
                    piece[-1] = '*'
                    termination = 1
                    piece.append((bset.get_block(block.no+1).t_start - 
                        (block.t_start+block.size)) * '_')
            
            
            else:
                # If this is the last block and the query does not map to the 
                # end of the target reference (usually because the query is 
                # only a coding region and is therefore not the same length),
                # pad the end to ensure the lengths match.
                #
                # Q >========|
                #    ////////
                # T >========XX|
                #
                if block.t_start + block.size < bset.t_size:
                        piece[-1] = '*'
                        piece.append((bset.t_size - 
                            (block.t_start + block.size)) * '_')

            pieces.append("".join(piece))
        final = "".join(pieces)

        if not args.quiet:
            #TODO Nicer output
            print "[READ]["+str(loopPos)+" of "+str(len(blocksets))+"]\t"+name.strip()
            print "\tFinal Sequence Length: "+str(len(final))
            print "======================================================================="

        if args.textwrap:
            output.write(textwrap.fill(final, 80))
        else:
            output.write(final)

        if not args.merged:
            output.close()
        else:
            output.write("\n")

        if not args.pslx:
            fasta_lib.close()
    blat_output.close()
    print "[DONE]\tAlignment Complete"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse .psl file output, "+
        "creating a new reference aligned FASTA file for each file processed"+
        " with blat.")
    parser.add_argument('blat_out', metavar="blat_out", type=str, 
        help=".psl output from BLAT run")
    parser.add_argument('-p', '--FASTA_path', metavar="", type=str,
        help="Path to extracted FASTA files directory processed by blat [default extracted/]",
        default="extracted/")
    parser.add_argument('-l', '--FASTA_list', metavar="", type=str,
        help="New line delimited list of FASTA file names processed by BLAT [default extracted/files.txt]",
        default="extracted/files.txt")
    parser.add_argument('-a', '--aligned', metavar="", type=str,
        help="Desired output directory for aligned FASTA files [default aligned/]",
        default="aligned/")
    parser.add_argument('-t', '--textwrap', action="store_true",
        help=("Output FASTA files are wrapped 80 characters a line, note that "
            "this can impact performance if there are a large "
            "number of sequences to process"))
    parser.add_argument('-q', '--quiet', action="store_true",
        help="Suppress read status output during FASTA processing")
    parser.add_argument('-x', '--pslx', action="store_true",
        help=("Read query sequence information for alignment from the pslx "
            "output instead of a FASTA file dir. Arguments for flags -p and "
            "-l will be ignored"))
    parser.add_argument('--library', metavar="", type=str,
        help=("Read from library FASTA instead of folder. Arguments for -p "
            "and -l will be ignored."))
    parser.add_argument('-m', '--merged', action="store_true",
        help=("Output all aligned sequences to the same file"))
    parser.add_argument('--multi', action="store_true",
        help=("MUST be set if you are providing multiple targets!"))#TODO? 
    parser.add_argument('--type', metavar="", type=str, default="psl/pslx",
        help=("Input file type, only does psl\pslx anyway [default psl/pslx]"))
    
    args = parser.parse_args()
    execute(args)
