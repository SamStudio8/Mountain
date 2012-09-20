#Mountain
A work in progress collection of simple but possibly useful genetics tools.

##Tools
###Aligner
Parse the output from an alignment tool (currently only blat) and create reference aligned FASTA files for each sample.<br />
Provides a library to read from the PSL/PSLX file format.

###Viewer
Curses based genome viewer, can load a folder or library of FASTA sequences and allows for traversal across the samples.
<br />Known SNP locations or points of interest can be highlighted given a file of new line delimited base locations.

###Leveller
Individually compares a folder or library of FASTA sequences base-by-base to a given reference and attempts to locate regions that could be eligible for assay for a given list of SNP locations.

##Contributors
* Sam Nicholls &lt;sn8@sanger.ac.uk&gt;

##License
Mountain is distributed under the GNU Lesser General Public License (LGPL 3).
See LICENSE.txt

##Copyright
Copyright &copy; 2012 Genome Research Ltd.
