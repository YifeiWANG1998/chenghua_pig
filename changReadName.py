import sys
import pysam

suffix = sys.argv[1]
infile = pysam.AlignmentFile('-', 'r')
outfile = pysam.AlignmentFile(sys.argv[2], 'wb', template=infile)
for idx, read in enumerate(infile):
	read.query_name = read.query_name + suffix
	outfile.write(read)
infile.close()
outfile.close()
