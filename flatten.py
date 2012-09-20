#!/usr/bin/env python
"""Individually compares a folder or library of FASTA sequences base-by-base
against a given reference and attempts to locate regions eligible for assay
for a given list of SNP locations."""

__author__ = "Sam Nicholls <sn8@sanger.ac.uk>"
__copyright__ = "Copyright (c) 2012 Genome Research Ltd."
__version__ = 0.850
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

from bs4 import BeautifulSoup as Soup
import argparse
import re
import os
import urllib2
import sys
import math

from Mountain.IO import fastareaders
from Mountain.Const import bases

class Flattener:
    def __init__(self):
        self.PHYLOGENY_URL = "http://www.mtdnacommunity.org/downloads/mtDNAPhylogeny.xml"
        self.PARSE_VARIANT = re.compile("([GCTAgcta])(\d+)([GCTAgcta])")

    def execute(self, args):
        self.FLANK_SIZE = args.flank
        self.MIN_FLANK = args.min_flank
        self.SNP_FLANK = args.oligo_flank
        self.MIN_OLIGO = args.min_oligo

        use_soup = 0 # Whether to read variant information using BeautifulSoup
        variants = []
        seen_variants = []
        if args.phylogeny:
            use_soup = 1
            print "[READ]\tReading mtDNA Phylogeny Tree from '"+args.phylogeny+"'"
            phylogeny_xml = open(args.phylogeny, 'r')
            print "[OKAY]\tmtDNA Phylogeny Tree Imported Successfully"
        elif args.list_variants:
            print "[READ]\tReading variant list from '"+args.list_variants+"'"
            variant_list = open(args.list_variants, 'r')
            for v in variant_list:
                if v not in seen_variants:
                    seen_variants.append(v.strip())
            print "[OKAY]\t"+str(len(seen_variants))+" Variant(s) Read"
        else:
            use_soup = 1
            phylogeny_xml = self.get_phylogeny(self.PHYLOGENY_URL)

        if use_soup:
            #TODO Should abstract file reading
            seen_variants = self.read_soup(phylogeny_xml)

        # Initalise dict for each detected variant
        variants = filter(None, map(self.parse_variant, seen_variants))
        
        # Read each base of the reference (target) sequence and initialise an
        # entry in a dictionary for each position.
        seqDict = self.init_refseq_dict(args.ref)

        fileChanges = {}
        currentPos = 0
        excluded = 0
        # Open each FASTA file in extracted directory.
        if os.path.isdir(args.aligned):
            use_library = False
            fastaList = sorted(os.listdir(args.aligned))
        else:
            use_library = True
            library = {}
            fa_library = fastareaders.FASTA_Library(args.aligned)
            for book in fa_library:
                library[book.readname] = book.seq
            fa_library.close()
            fastaList = sorted(library.keys())

        fileCount = len(fastaList)
        for fasta in fastaList:
            currentPos += 1
            fileChanges[fasta] = 0

            stro = "[READ]["+str(currentPos)+" of "+str(fileCount)+"]\t"+fasta+"\r"
            sys.stdout.write(stro)
            sys.stdout.flush()
            sys.stdout.write((len(stro)*" ")+"                \r")

            if use_library:
                seq = library[fasta]
            else:
                handler = fastareaders.FASTA_Library(args.aligned+fasta)
                fa = handler.get_next()
                seq = fa.seq

            # Exclude FASTA files that are not from the M or N phylogeny branches
            if args.onlyMN:
                try:
                    if (seq[10399] == "T" or seq[9539] == "T"):
                        pass
                    else:
                        excluded += 1
                        continue
                except IndexError:
                    excluded += 1
                    continue

            for i, n in enumerate(seq):
                n = n.upper()
                if n != "_" and n != "*" and n != "N":
                    if n in seqDict[i]['base_counter']:
                        seqDict[i]['base_counter'][n] += 1
                    else:
                        # Non standard base, probably IUPAC
                        seqDict[i]['base_counter'][n] = 1
                elif n == "*" or n == "N":
                    # Ignore '_' completely.
                    seqDict[i]['base_counter']["N*"] += 1
                    
                # Ignore 'N's in query sequence
                if args.N:
                    if n == "N":
                        continue

                # PSL Alignment tool now uses "_" and marks termination 
                # points with * making the following block redundant
                # Ignore missing bases in query sequence
                #if args.nodel:
                #    if n == "-":
                #        continue

                if (seqDict[i]['base'] != 'N' 
                        and n != seqDict[i]['base']
                        and n != '_'):

                    # Try to reconcile non standard bases to a IUPAC code
                    if args.iupac:
                        if n in bases.IUPAC:
                            if seqDict[i]['base'] in bases.IUPAC[n]:
                                continue

                    #TODO Flag to generate exclusion list.
                    #TODO Parameter to accept optional exclusion list.
                    seqDict[i]['count'] += 1
                    fileChanges[fasta] += 1
                    if seqDict[i]['count'] > args.tolerance:
                        seqDict[i]['base'] = 'N'

            if not use_library:
                handler.close()
        
        # Construct Sequence
        converted = 0
        processedSeq = []
        for k in sorted(seqDict.keys()):
            if args.iupac_out:
                current_base_failure_threshold = math.ceil(
                    math.ceil(100 * (float(args.tolerance)/float(fileCount-excluded))) *
                    (float(sum(seqDict[k]['base_counter'].values())) / float(100)))

                # "Upbase" the standard base count for any detected IUPAC codes
                #TODO Double check this - (1k data didn't need it)
                for intersected_code in (set(bases.IUPAC) & set(seqDict[k]['base_counter'])):
                    for standard_base in bases.IUPAC[intersected_code]:
                        seqDict[k]['base_counter'][standard_base] += 1

                # "Upbase" the standard base count for any N's or *'s detected
                if "N*" in seqDict[k]['base_counter']:
                    for standard_base in bases.STANDARD_BASES:
                        seqDict[k]['base_counter'][standard_base] += seqDict[k]['base_counter']["N*"]

                current_base_threshold_bases = dict(filter(
                    lambda (key, val) : val > current_base_failure_threshold
                    and key != "N*",
                    [(key, val) for key, val in seqDict[k]['base_counter'].items()]))

                current_iupac_flag = 0
                for iupac in bases.IUPAC:
                    if sorted(current_base_threshold_bases) == sorted(bases.IUPAC[iupac]):
                        if seqDict[k]['base'] != "N":
                            print k
                            print seqDict[k]['base']
                            print seqDict[k]['base_counter']
                        converted += 1
                        current_iupac_flag = 1
                        processedSeq.append(iupac)

                if not current_iupac_flag:
                    processedSeq.append(seqDict[k]['base'])
            else:
                processedSeq.append(seqDict[k]['base'])
        
        if args.blocks:
            blockout_block = 0
            blockout_line = 0
            output = open(args.out+'.'+str(blockout_block), "w+")
        else:
            output = open(args.out, "w+")

        possibleProbesCount = 0
        for v in variants:
            fin = self.make_probe(v, processedSeq)
            probe_flanks = self.validate_probe(fin)
            if probe_flanks:
                possibleProbesCount += 1
                if args.primer3:
                    output.write("SEQUENCE_ID="+v['variant']+"\n")
                    output.write("SEQUENCE_TEMPLATE="+_probe_flanks[0]+v['x']+probe_flanks[2]+"\n")
                    output.write("SEQUENCE_TARGET="+str(len(probe_flanks[0])+1)+",1"+"\n")
                    output.write("=\n")
                else:
                    curv = re.sub("!", "bang", v['variant'])
                    output.write(">"+curv+"\n")
                    output.write(fin+"\n")

                    if args.blocks:
                        blockout_line += 1
                        if blockout_line == args.blocks:
                            output.close()
                            blockout_block += 1
                            blockout_line = 0
                            output = open(args.out+'.'+str(blockout_block), "w+")
                    
    #    counterOut = open("counter", "w+")
    #    for i in seqDict:
    #        counterOut.write(str(i+1)+"\t"+str(seqDict[i]['count'])+"\n")
    #    counterOut.close()

    #    fastaWarningCount = 0
    #    fileOut = open("fileReport", "w+")
    #    for f in sorted(fileChanges.iteritems(), key=lambda x: x[1], reverse=True):
    #        fileOut.write(str(f)+"\n")
    #        if f[1] >= 100:
    #            fastaWarningCount += 1
    #    fileOut.close()
        
    #    if fastaWarningCount > 0:
    #        print "[WARN]\t"+str(fastaWarningCount)+" FASTA files with more than 100 bases different from reference"

        if args.onlyMN:
            print "[NOTE]\t"+str(excluded)+" FASTA Files Excluded"

        if args.iupac_out:
            print "[NOTE]\t"+str(converted)+" Bases Converted to IUPAC Codes"

        print "[DONE]\t"+str(possibleProbesCount)+" Eligible Sequences Detected"
        output.close()

    def get_phylogeny(self, location):
        """Download and return XML mtDNA Phylogenetic Tree"""

        print "[HTTP]\tFetching "+location
        url = urllib2.urlopen(location)
        print "[READ]\tReading mtDNA Phylogeny Tree"
        phylogeny_xml = url.read()
        url.close()
        print "[OKAY]\tmtDNA Phylogeny Tree Imported Successfully"
        return phylogeny_xml

    def read_soup(self, phylogeny):
        """Parse XML mtDNA Phylogenetic Tree with BeautifulSoup and return a list
        of known mtDNA variants"""

        #TODO Ignore bangs and parens
        seen_variants = []
        print "[SOUP]\tParsing Phylogeny Tree XML"
        soup = Soup(phylogeny, "xml")
        for node in soup.find_all(HG=True):
            for v in re.split(",", node["HG"]):
                if v.strip() not in seen_variants:
                    seen_variants.append(v.strip())
        print "[OKAY]\t"+str(len(seen_variants))+" Variant(s) Read"
        return seen_variants

    def init_refseq_dict(self, ref_path):
        #TODO seqDict[i][isSNP] = [1|0]
        """Load and traverse the reference sequence base-by-base, creating a dict
        entry with the 0-index base location key with dict as value containing:
            the reference base itself, 
            a "mutation" threshold counter and 
            a dict that keeps track of all bases seen at this location"""

        seqDict = {}
        ref_fasta = fastareaders.FASTA_Library(ref_path).get_next()
        for i, base in enumerate(ref_fasta):
            base = base.upper()

            # Warn for non standard bases as query sequences will not match to
            # an IUPAC or 'N' base in the reference.
            if re.match("[^GCTA]", base):
                print ("[WARN]\tNon GCTA Nucleotide in Reference at 0-index: "
                    +str(i)+" ("+base+")")

            seqDict[i] = {}
            seqDict[i]['base'] = base
            seqDict[i]['count'] = 0
            seqDict[i]['base_counter'] = {'G':0, 'T':0, 'C':0, 'A':0, 'N*':0}

        return seqDict

    def parse_variant(self, v):
        """Parse the name of an mtDNA variant and return a dictionary containing
        the variant name, the x and y allele and the 1-indexed base position"""

        groups = self.PARSE_VARIANT.search(v)
        if groups:
            return {
                'variant': v,
                'x': groups.group(1), # x allele
                'y': groups.group(3), # y allele
                'location': int(groups.group(2)) # 1-indexed base position
            }
        else:
            return None

    def validate_probe(self, seq):
        """Check whether a given subsequence could be an eligible probe"""
        #
        #  Left Flank        <<|X|>>       Right Flank
        #  --------------------| |--------------------
        #  >====================X====================>
        #                       ^SNP Site
        #
        #  The "flanks" to the left and right of the SNP location, which 
        #  extend all the way to the beginning and end of the subsequence 
        #  must *both* contain at least one straight run of standard bases 
        #  of length args.min_flank anywhere in the flanking region.
        #
        #
        #
        #                        Internal Oligo
        #                    ---------------------
        #  >================|===X=================|==>
        #  >=========|==========X==========|=========>
        #             ---------------------
        #                Internal Oligo
        #
        #  For the "internal oligo" to exist, it must consist of at least 
        #  one straight run of standard bases of length args.min_oligo, 
        #  and must pass through the SNP site.
        #  Note that it need not be symmetrical to the SNP site but there
        #  must be at least one base either side of the SNP in the sequence.
        #
        #
        #
        #  >==================--X--==================>
        #                     ^^ ^^
        #
        #  Finally, there must be no non-standard bases either side of the
        #  SNP site for a run of at least args.oligo_flank bases.
        #
        SPLIT_PROBE = re.compile("(\[\w/\w\])")
        OLIGO = re.compile("([GCTAgcta]{1,})(\[\w/\w\])([GCTAgcta]{1,})")
        FLANK = re.compile("[GCTAgcta]{"+str(self.MIN_FLANK)+",}")
        CENTER_CHECK = re.compile("[GCTAgcta]{"+str(self.SNP_FLANK)+"}"+
            "\[\w/\w\][GCTAgcta]{"+str(self.SNP_FLANK)+"}")

        matches = OLIGO.search(seq)
        if not matches:
            return False
        if (len(matches.group(1)) + len(matches.group(3)) + 1) >= self.MIN_OLIGO:
            flanks = SPLIT_PROBE.split(seq)
            left = FLANK.search(flanks[0])
            right = FLANK.search(flanks[2])
            center = CENTER_CHECK.search(seq)
            if left and right and center:
                return flanks
        else:
            return False

    def make_probe(self, v, processedSeq):
        """Return the probe region string for a given variant with flanks of 
        the given size"""
        if v['location'] < self.FLANK_SIZE:
            start = 0
            end = v['location']+self.FLANK_SIZE

            region = processedSeq[0 : v['location']+self.FLANK_SIZE]
            offset = self.FLANK_SIZE-v['location']+1
        else:
            start = v['location']-1-self.FLANK_SIZE
            end = v['location']+self.FLANK_SIZE

            region = processedSeq[v['location']-1-self.FLANK_SIZE : v['location']+self.FLANK_SIZE]
            offset = 0
        region[self.FLANK_SIZE-offset] = "["+v['x']+"/"+v['y']+"]"
        return "".join(region)

if __name__ == '__main__':
    #TODO Add description and improve help comments
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    parser.add_argument('ref', metavar="reference", type=str,
        help=("The reference FASTA used for alignment"))
    parser.add_argument('out', metavar="out", type=str,
        help="Output location")
    group.add_argument('-p', '--phylogeny', metavar="", type=str,
        help=("For all variants: Latest version of XML mtDNA Phylogenic Tree "
            "[default: Automatic download from mtdnaCommunity.org]"))
    group.add_argument('-l', '--list_variants', metavar="", type=str,
        help=("For particular variants: Path to new line delimited file of "
            "mtDNA variants for processing. If this option is not selected, "
            "the script will default to downloading the phylogenic tree."))
    parser.add_argument('--N', action="store_true",
        help=("'N' base in query sequences will match any nucleotide in "
            "reference"))
    parser.add_argument('--iupac', action="store_true",
        help=("Attempt to match non-standard bases found in query sequences "
            "to an IUPAC code that contains the reference base"))
#    parser.add_argument('--nodel', action="store_true",
#        help=("Deletions found in the query sequence are still matched to "
#            "the target reference"))
    parser.add_argument('--onlyMN', action="store_true",
        help=("Query sequences not in the M or N haplotype are ignored"))
    parser.add_argument('--min_oligo', metavar="", type=int,
        help=("Minimum sequence size containing 0 N bases that must run "
            "through the SNP [default: 13]"),
        default=13)
    parser.add_argument('--min_flank', metavar="", type=int,
        help=("Minimum subsequence size containing 0 N bases that must flank "
            "the SNP somewhere in the sequence [default: 9]"),
        default=9)
    parser.add_argument('-f', '--flank', metavar="", type=int,
        help=("Length of flanking sequences either side of the SNP base to"
            "process [default: 400]"),
        default=400)
    parser.add_argument('-t', '--tolerance', metavar="", type=int,
        help=("Mutation threshold for each base, ie. the number of query "
            "sequences that are allowed to differ from the reference at "
            "each base [default: 0]"),
        default=0)
    parser.add_argument('-p3', '--primer3', action="store_true",
        help=("Format output for Primer3 [default: FASTA]"))
    parser.add_argument("--iupac_out", action="store_true",
        help=("Output an IUPAC code instead of N where possible"))
    parser.add_argument("--oligo_flank", type=int, default=2)
    parser.add_argument("-a", "--aligned", metavar="",
        help=("FASTA folder or Library"))
    parser.add_argument('--blocks', type=int, metavar=""
        help=("Output results in blocks, with n sequences per block"))
    #TODO Flag to output SNPs as N regardless of what was read
    Flattener().execute(parser.parse_args())
