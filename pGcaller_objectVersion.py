#! /bin/env python
import time
startime = time.time()
import sys
import argparse
import subprocess
from cStringIO import StringIO
import math
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

############################################
# This script requires python 2.7 or later #
############################################
#20110711 Erik Borgstroem


#commandline argument parsing
argparser = argparse.ArgumentParser(	description='Script for calling length of polyG regions in reads extracted from a ".bam" file. The scripts require that samtools are installed (at uppmax run "module load samtools").',
										formatter_class=argparse.RawTextHelpFormatter,
									)
argparser.add_argument(	'-b',	dest='bamfile',					metavar='FILE',	type=file,				required=True,					help='Input ".bam" file to extract reads from. Must use full path!')
argparser.add_argument(	'-p',	dest='pgfile',					metavar='FILE',	type=file,				required=True,					help='Fasta file containing polyG regions found in hg19/grch37.')
argparser.add_argument(	'-o',	dest='outfile',					metavar='str',	type=str,				required=True,					help='Output filename prefix, file were info about each polyG will be printed.')
argparser.add_argument(	'-m',	dest='match',					metavar='N',	type=int,				required=False,	default=1,		help='Number of bases that need to match at each side of polyG region eg, -m 3 means N(G*n)N (defaults to minimum value = 1).')
argparser.add_argument(	'-q',	dest='min_mapq',				metavar='N',	type=int,				required=False,	default=0,		help='Minimum mapping quality requiered to consider read.')
argparser.add_argument(	'-l',	dest='load_max',				metavar='N',	type=int,				required=False,	default=0,		help='Maximum number of pg reagions to load, for testing & debugging purposes (0 = inifint, default).')
argparser.add_argument(	'-pci',	dest='print_calling_info',						action='store_true',	required=False,	default=False,	help='Print information during length calling to output (default= Flase).')
argparser.add_argument(	'-pma',	dest='print_multi_alignemnt',					action='store_true',	required=False,	default=False,	help='Print multi alignemnt of reads covering each polyG (default= Flase).')
argparser.add_argument(	'-nbb',	dest='no_bold_bases',							action='store_true',	required=False,	default=False,	help='Do not show mismatch bases in multialignment with bold font (default= false).')
argparser.add_argument(	'-psc',	dest='print_samtools_command',					action='store_true',	required=False,	default=False,	help='Print the samtools command used for read extraction from bamfile, for debugging purposes (default= false).')
argparser.add_argument(	'-rrc',	dest='required_read_count',		metavar='N',	type=int,				required=False,	default=4,		help='The read count requiered for calling of polyG length (default = 4).')
argparser.add_argument(	'-ps',	dest='percentage_support',			metavar='%',	type=int,				required=False,	default=40,		help='Percentage support needed to call an allel (default = 40, must be >1/3 or three allels can be called).')
# percentage support requiered for heterozygote call
input = argparser.parse_args(sys.argv[1:])
assert input.match != 0, "atleast one base has to match on each side of polyG to be able do do analysis, ie -m cannot be zero."
#convert and prepare input 
input.percentage_support = input.percentage_support/100.00
assert input.percentage_support > 1.0/3.0, 'percentage_support (-ps) has to be larger than 1/3'
# strings & settings for bold output
b_start = "\033[1m"
b_end = "\033[0;0m"
	

def main():
#--------------------------MAIN PROGRAM-----------------------------
	pass
#--------------------------MAIN PROGRAM END-------------------------


#--------------------- Functions // Subroutines --------------------
class reference():
	"""Object that holds reference info."""
	def __init__(self, length,seq):
		self.length	= length
		self.seq	= seq




class PolyG():
	"""The polyG class holds info about a polyG region and it's reads aswell as functions for calling pg lengths and statistics."""

	def __init__(self,pgid,chrom,pos,length,seq):
		#information from fasta database:
		self.pgid	= pgid
		self.chrom	= chrom
		self.pos	= pos
		self.reference = reference(length,seq);
		#Flags to keep track of actions performed:
		self.alignment_str_flag = False;
		self.reads_filtered = False;
		self.reads_loaded = False;
		self.call_flag = 'NotStarted'
		#information to be filled:
		self.output = ''
		self.read_count = 0
		self.filtered_read_count = 0
		self.pre_filter_read_count = None




	def getReads(self):
		"""Function that extract reaads from bamfile and stores them as 'read' objects in self.reads."""

		class read():
			"""Object that holds an aligned sequence read."""
			def __init__(self, line):
				line = line.split('\t')
				try:
					self.qname	= line[0] #Query (pair) NAME
					self.flag	= int(line[1]) #bitwise FLAG
					self.rname	= line[2] #Reference sequence NAME
					self.pos	= int(line[3]) #1-based leftmost POSition/coordinate of clipped sequence
					self.mapq	= int(line[4]) #MAPping Quality (Phred-scaled)
					cigar_pattern = re.compile('[0-9]+[MIDNSHP]')
					self.cigar	= cigar_pattern.findall(line[5]) #extended CIGAR string
					self.mrnm	= line[6] #Mate Reference sequence NaMe ("=" if same as RNAME)
					self.mpos	= line[7] #1-based Mate POSistion
					self.isize	= line[8] #Inferred insert SIZE
					self.seq	= line[9] #query SEQuence on the same strand as the reference
					self.qual	= line[10] #query QUALity (ASCII-33 gives the Phred base quality)
					self.opt	= line[11:] #variable OPTional fields in the format TAG:VTYPE:VALUE
					self.length	= len( self.seq)
					assert self.rname	!= '*'
					assert self.pos		!= 0
				except IndexError:
					sys.stderr.write( "IndexError!!"+'\n')
					sys.stderr.write( str(line)+'\n')
					sys.exit()
		#End of read class

		# ---- getReads() function work starting here ----

		#Store command for debugging purposes
		self.get_read_cmd = ' '.join(['samtools', 'view', input.bamfile.name, self.chrom+':'+str(self.pos)+'-'+str( self.pos + self.reference.length -1)])

		#startsubprocess that extracts samfile viewing the specified region
		samtools = subprocess.Popen(['samtools', 'view', input.bamfile.name, self.chrom+':'+str(self.pos)+'-'+str( self.pos + self.reference.length -1)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		samfile, errdata = samtools.communicate()
		if samtools.returncode != 0:
			print 'Samtools view Error code', samtools.returncode, errdata
			sys.exit()
		samfile = StringIO(samfile)

		#go through each line in samfile and add read to array
		self.reads = []
		for line in samfile:
			self.reads.append(read(line))
		self.read_count = len(self.reads)

		#set flag to True
		self.reads_loaded = True
	#end of getReads()




	def filter_reads(self):
		"""Removes reads that do not cover the polyG regionor has to low mapping quality."""
		if not self.reads: return # if reads are empty do nothing
		temporary = []
		for read in self.reads: #check read positions
			trimmed_length = read.length
			for block in read.cigar: #check cigar string and decrease trimmed read length
				last_block=False
				first_block=False
				if block == read.cigar[-1]:last_block=True
				if block == read.cigar[0]:first_block=True
				bases = int(block[:-1]) # block length
				type = block[-1:]		# block type
				if type == 'S':
					if last_block or first_block:
						trimmed_length -= bases
					else: sys.stderr.write('Error: S has to be first or last block in CIGAR')
			if read.pos > self.pos-input.match: continue #Skips if read starts within pg region or edge bases
			if read.pos + trimmed_length < self.pos + self.reference.length + input.match: continue #Skips if read ends within pg region or edge bases
			if read.mapq < input.min_mapq: continue #Mapping quality below threshold => skip read
			temporary.append(read) #append to temp filtered reads array if all test passed
		self.reads = temporary #set self.reads to filtered reads array
		self.filtered_read_count = len(self.reads)
		self.reads_filtered = True
	#end of filter_reads()




	def get_pre_filter_read_count(self):
		"""Count the number of reads that cover the full polyG region."""
		self.pre_filter_read_count = 0
		if not self.reads:return # if reads are empty do nothing
		for read in self.reads:
			if read.pos > self.pos: continue #Skips if read starts within pg region
			if read.pos + read.length < self.pos + self.reference.length: continue #Skips if read ends within pg region
			self.pre_filter_read_count += 1
	#end of filter_reads()




	def make_alignemnt_strings(self):
		"""Create strings used during printing of alignemnt view"""
#		self.alignment_str_flag = '0:started'
		for read in self.reads:
			temporary = ''
			current_pos = 0
			for part in read.cigar:
				bases = int(part[:-1])
				type = part[-1:]
				if type == 'M': #Match or Missmatch view part in alignement
					temporary += read.seq[current_pos:current_pos+bases]
					current_pos += bases
				elif type == 'I': #Insertion not viewed in alignment
					current_pos += bases
				elif type == 'D': #Deletion viewed as '-' in alignment
					temporary += ''.join(['-' for i in range(bases)])
				elif type == 'S': #"remove" soft trimmed bases from alignment string
					current_pos += bases
				#currently unsupported types in cigar string
				elif type == 'N': sys.stderr.write( 'Error: '+type+' unsupported in CIGARstring!!!\n');sys.exit() #Skipped
				elif type == 'P': sys.stderr.write( 'Error: '+type+' unsupported in CIGARstring!!!\n');sys.exit() #Padding
				elif type == 'H': sys.stderr.write( 'Error: '+type+' unsupported in CIGARstring!!!\n');sys.exit() #Hard clipping
			read.alignment_str = temporary
		self.alignment_str_flag = True
	# end of make_alignemnt_strings()




	def print_alignment(self):
		"""Prints an alignment of polyG and reads, note that insertions are not displayed other than by the cigar string"""
		self.make_alignemnt_strings() # make alignment strings
		self.output += ( self.reference.seq+' <-ref\tdir\tcigar\n')

		for read in self.reads:
			print_start	= 0
			print_end	= len(read.alignment_str)
			ref_start	= 0
			ref_end		= len(self.reference.seq)
			ref_overhang= 0

			#calculate relative position to reference sequence and print spaces (-) if needed before read start
			if read.pos > self.pos-50:
				if print_end + (read.pos-(self.pos-50)) > len(self.reference.seq): print_end -= print_end + (read.pos-(self.pos-50)) - len(self.reference.seq) 
				ref_start = (read.pos-(self.pos-50))
				self.output += (''.join(['-' for i in range(read.pos-(self.pos-50))]))
			if read.pos < self.pos-50:
				if print_end - (self.pos-50-read.pos) > len(self.reference.seq): print_end -= print_end - (self.pos-50-read.pos) - len(self.reference.seq) 
				print_start = (self.pos-50-read.pos)
			if ref_start + (print_end - print_start) < ref_end:
				ref_overhang = (ref_end - (ref_start + (print_end - print_start)))
				ref_end = ref_end - ref_overhang

			#find mismatches and print in bold font
			for base in zip(read.alignment_str[print_start:print_end],self.reference.seq[ref_start:ref_end]):
				if base[0] == base[1]: self.output += (base[0])
				else:
					if input.no_bold_bases: mmbase_string = base[0];
					else: mmbase_string = b_start+base[0]+b_end
					self.output += (mmbase_string)
			
			#print spaces after read if needed (and other info)
			self.output += (''.join(['-' for i in range(ref_overhang)]))
			self.output += ('\t'+explainflag(read.flag)[4])
			self.output += ('\t'+''.join(read.cigar))
#			self.output += ('\t'+str(read.pos))
			self.output += ('\n')
			assert len(read.alignment_str[print_start:print_end]) == len(self.reference.seq[ref_start:ref_end])
		self.output += ('\n') 
		self.cleanup(3) #remove alignemnt strings
	#end of print_alignment()




	def call_pglengths(self):
		"""Function that calls polyG lengths in reads"""
		#Set initial values and check that reads are present
		self.sucessfully_called_reads = 0
		first_read = True
		if not self.reads: sys.stderr.write('Warning: '+self.pd_id+' has no reads and no calling has been done');self.call_flag = '1:no reads present';return
		else: self.call_flag = '2:readcalling started'

		for read in self.reads:					# for ech read present do the following:
			read.call_flag = '0:OK'				# set flag calling OK
			read.pg_pos = self.pos - read.pos	# calculate position of polyGregion in the read
													# note that this is relative distance (only number of bases between readpos and pg startpos) ie. first base has position zero, in cigar first pos is one!!!!
			#set default read flags:
			D_in_pg		= False
			I_in_pg		= False
			I_before_pg	= False
			read.MM_in_pg	= False
			l_edge_D	= False
			r_edge_D	= False
			l_edge_I	= False
			r_edge_I	= False
			edge_MM		= False
			ref_position = read.pos #our current position in the reference
			read.called_length = self.reference.length

			#Go through each block in the cigar string and look for mm, insertions or deletions in polyGregion
			for block in read.cigar:
				last_block=False
				first_block=False
				if block == read.cigar[-1]:last_block=True
				if block == read.cigar[0]:first_block=True
				bases = int(block[:-1]) # block length
				type = block[-1:]		# block type
				check_end = False

				#calculate block start position in reference relative to polyG
				start = False
				if ref_position		<=	self.pos-2: start = 'before'					#block starting >2 base before polyG startpos
				elif ref_position	==	self.pos-1:	start = 'l_edge_m1'					#block starting 2 bases before polyG startpos
				elif ref_position	==	self.pos: start = 'l_edge'						#block starting one base before polyG startpos
				elif ref_position	<	self.pos+self.reference.length:	start = 'within'#block starting from on first base of polyG to last base of polyG
				elif ref_position	==	self.pos+self.reference.length: start = 'r_edge'#block one base after polyG
				elif ref_position	>	self.pos+self.reference.length:	start = 'after'	#block starting >1 base after polyG end pos
				else: raise ValueError, 'position is not possible'

				#check type of block and adjust position of polyG in read (and increment reference positions for next block), following pg will be used for explanation CAGTCCCCCTGAC
				if type == 'M': #Match or Missmatch, same number of bases ie. increase reference positions with "bases", no shift in polyG postition
					ref_position += bases

				elif type == 'I': #Insertion, 'bases' base more in read than reference ie. shift the polyG position forward (increase read.pg_pos), same reference position
					if start == 'before' or  start == 'l_edge_m1':	# insertion starts before pg position ie just shift pg pos (eg [N*n]CAGTCCCCCTGAC), if close to pg eg. CAG[N*n]TCCCCCTGAC this will lead to edge mismatch
						read.pg_pos += bases
						I_before_pg = True
					elif start == 'l_edge':	# the insertioin starts after the first base before pg ie edge bases will match (CAGT[N*n]CCCCCTGAC), is insertion pg region or not? check sequence? FIX!?
						read.called_length += bases
						I_in_pg = True
						l_edge_I = True
					elif start == 'within':	# the insertion starts after the polyG start position and within pg region but not after the last pg position (CAGTC[N*n]C[N*n]C[N*n]C[N*n]CTGAC)
						read.called_length += bases
						I_in_pg = True
					elif start == 'r_edge':	# the insertion starts after the last base in the pg ie edge bases will match (CAGTCCCCC[N*n]TGAC), is insertion pg region or not? check sequence? FIX!?
						read.called_length += bases
						I_in_pg = True
						r_edge_I = True
					elif start == 'after':	# the insertion starts after the pg ie either will the edge bases not match or the calling will be uneffected (CAGTCCCCCT[N*n]GAC to CAGTCCCCCTGAC[N*n])
						pass
					else: raise ValueError

				elif type == 'D': #Deletion, 'bases' base less in read than in reference  ie. increase ref position, shift polyG position backward in read (decrease read.pg_pos)
					ref_position += bases
					if start == 'before':	# deleteion starts before pg position ie just shift pg pos (eg [-*n]CAGTCCCCCTGAC)
						read.pg_pos -= bases# if del is close to pg eg. CA-TCCCCCTGAC this will lead to edge mismatch
						check_end = True	# when longer than 1bases check if end pos is base before pg, in pg or after polyG (CA--CCCCCTGAC or -------CCTGAC or CA---CCCCTGAC to CA---CCCCTGAC or CA--------GAC) calling might be unvalid as >=1 base on each side have to match, flag for this 
					elif start == 'l_edge_m1':# deleteion starts the base before the pg, ie edge bases will not match (CAG-CCCCCTGAC or CAG--CCCCTGAC or CAG-------GAC)
						read.pg_pos -= bases# no need check end position as calling is unvalid in any case as >=1 base on each side have to match
						l_edge_D = True
					elif start == 'l_edge':	# deleteion starts on the first base of the pg (CAGT-CCCCTGAC), check if end position is in or after polyG (CAGT---CCTGAC=valid or CAGT------GAC=unvalid)
						read.called_length -= bases
						D_in_pg = True
						check_end = True
					elif start == 'within':	# the deletion starts from 2nd base of pG to last base of pG (CAGTC-CCCTGAC to CAGTCCCC-TGAC), check if end position is in or after polyG (CAGTC--CCTGAC=valid or CAGTC-----GAC=unvalid)
						read.called_length -= bases
						D_in_pg = True
						check_end = True
					elif start == 'r_edge':	# the deletion starts on the first base after the pG (CAGTCCCCC-GAC or CAGTCCCCC--AC and so on...), this is always unvalid as >=1 base on each side have to match
						r_edge_D = True
					elif start == 'after':	# the Deletion starts after the pg ie either will the edge bases not match or the calling will be uneffected (CAGTCCCCCT-AC to CAGTCCCCCTGAC[-*n]), will most likely be caught by the edge mismatch check and will otherwiae not affect the calling
						pass
					else: raise ValueError

				elif type == 'S':
					# pg34210
					# S blocks are always is in 3 or 5 terminal of read
					# skip this sequence in all calculations ie pretend that read is shorter
					# Read mapping position start with first base in Cigar that has M,IorD ie the S parts will not affect the reference positions but only the position of the polyG in the read
					# therefore if the S part overlaps with edge bases or polyG the read will be filtered away in the filter_reads() function
					# tough we still need to compensate the read pg position for all first_block trimming
					if first_block: # move pg in read position towards reference 3' terminal
						read.pg_pos += bases;
					elif last_block:
						pass
					else: sys.stderr.write('Error: '+type+' Has to be first or last in CIGAR!!!'+'\n');sys.exit() #Soft trimming

				#currently unsupported types in cigar string
				elif type == 'N': sys.stderr.write('Error: '+type+' unsupported in CIGARstring!!!'+'\n');sys.exit() #Skipped
				elif type == 'P': sys.stderr.write('Error: '+type+' unsupported in CIGARstring!!!'+'\n');sys.exit() #Padding
				elif type == 'H': sys.stderr.write('Error: '+type+' unsupported in CIGARstring!!!'+'\n');sys.exit() #Hard clipping

				if check_end: #Check end of deletion
					#calculate block end position in reference
					end = False
					if ref_position		<=	self.pos-2: end = 'before'						#block ending >2 base before polyG startpos
					elif ref_position	==	self.pos-1:	end = 'l_edge_m1'					#block ending 2 bases before polyG startpos
					elif ref_position	==	self.pos: end = 'l_edge'						#block ending one base before polyG startpos
					elif ref_position	<	self.pos+self.reference.length:	end = 'within'	#block ending from on first base of polyG to last base of polyG
					elif ref_position	==	self.pos+self.reference.length: end = 'r_edge'	#block one base after polyG
					elif ref_position	>=	self.pos+self.reference.length:	end = 'after'	#block ending >1 base after polyG end pos
					else: sys.stderr.write('refpos'+str(ref_position)+' pgpos'+str(self.pos)+' reflength'+str(self.reference.length)+'\n'); raise ValueError, 'end position is not possible'
					#check end position
					if start == 'before':
						if		end == 'before': pass
						elif	end == 'l_edge_m1': pass
						elif	end == 'l_edge' or end == 'within' or end == 'r_edge' or end == 'after': l_edge_D = True
						else: raise ValueError
					elif start == 'l_edge' or start == 'within':	
						if		end == 'before' or end == 'l_edge_m1' or end == 'l_edge' or end == 'within' or end == 'r_edge':	pass
						elif	end == 'after': r_edge_D = True
						else: raise ValueError
					else: raise ValueError

			#output info:
			if first_read and input.print_calling_info: self.output += ('id\t\t\t\tbefore:\tpolyG:\tafter:\tlength:\tflag:\ttags:\nreference sequence:       \t'+ self.reference.seq[50-input.match:50]+'\t'+self.reference.seq[50:50+self.reference.length]+'\t'+self.reference.seq[50+self.reference.length:50+self.reference.length+input.match]+'\t'+str(self.reference.length)+'\n');first_read=False
			if len(read.qname) < len('reference sequence:') and input.print_calling_info:self.output += ('\t'+''.join([' ' for i in range(len('reference sequence:')-len(read.qname))]))
			if input.print_calling_info: self.output += (read.qname+'\t')

			#look for mismatches in polyG and Edge bases (and print info)
			#the before polyG edge bases
			for base in zip(read.seq[read.pg_pos-input.match:read.pg_pos],self.reference.seq[50-input.match:50]):
				if base[0] == base[1]:
					if input.print_calling_info: self.output += (base[0])
				else:
					if input.print_calling_info:
						if input.no_bold_bases: mmbase_string = base[0];
						else: mmbase_string = b_start+base[0]+b_end
						self.output += (mmbase_string);
					if read.call_flag[0] != '5': read.call_flag = '1:Mismatch in edge bases';
			#within the polyG region
			if input.print_calling_info: self.output += ('\t')
			if read.called_length == self.reference.length:
				for base in zip(read.seq[read.pg_pos:read.pg_pos+self.reference.length],self.reference.seq[50:50+self.reference.length]):
					if base[0] == base[1]:
						if input.print_calling_info: self.output += (base[0])
					else:
						if input.print_calling_info:
							if input.no_bold_bases: mmbase_string = base[0];
							else: mmbase_string = b_start+base[0]+b_end
							self.output += (mmbase_string);
						read.MM_in_pg = True;
			else: 
				for base in zip(read.seq[read.pg_pos:read.pg_pos+read.called_length],''.join([self.reference.seq[50:51] for i in range(read.called_length)])):
					if base[0] == base[1]:
						if input.print_calling_info: self.output += (base[0])
					else: 
						if input.print_calling_info:
							if input.no_bold_bases: mmbase_string = base[0];
							else: mmbase_string = b_start+base[0]+b_end
							self.output += (mmbase_string);
						read.MM_in_pg = True				
			#the after polyG edge bases
			if input.print_calling_info: self.output += ('\t')
			for base in zip(read.seq[read.pg_pos+read.called_length:read.pg_pos+read.called_length+input.match],self.reference.seq[50+self.reference.length:50+self.reference.length+input.match]):
				if base[0] == base[1]:
					if input.print_calling_info: self.output += (base[0])
				else:
					if input.print_calling_info:
						if input.no_bold_bases: mmbase_string = base[0];
						else: mmbase_string = b_start+base[0]+b_end
						self.output += (mmbase_string);						
					if read.call_flag[0] != '5': read.call_flag = '2:Mismatch in edge bases';
			if input.print_calling_info: self.output += ('\t'+str(read.called_length))

			#set edge deletion flag
			if l_edge_D or r_edge_D: read.call_flag = '3:Deletion of edge base(s)'

			#print info
			if input.print_calling_info:
				self.output += ('\t'+read.call_flag[0])
				self.output += ('\t'+''.join(read.cigar))
				if D_in_pg: self.output += ('\tin_pG_Deletion')
				if I_in_pg: self.output += ('\tInsertion')
				if l_edge_I: self.output += ('\tl_edge_I')
				if r_edge_I: self.output += ('\tr_edge_I')
				if l_edge_D: self.output += ('\tl_edge_D')
				if r_edge_D: self.output += ('\tr_edge_D')
				if read.MM_in_pg: self.output += ('\tMisMatch')
				if read.call_flag[0] == '1': self.output += ('\tEdge_MM')
				elif read.call_flag[0] == '2': self.output += ('\tEdge_MM')
				elif read.call_flag[0] == '3': self.output += ('\tEdge_Deleted')
				elif read.call_flag[0] == '4': self.output += ('\tNoOverlap(shift)')
				elif read.call_flag[0] == '5': self.output += ('\tNoOverlap(Deleted)')
				self.output += ('\n')
			read.call_flag += ':finished'
			if read.call_flag[0] == '0': self.sucessfully_called_reads += 1
		if input.print_calling_info: self.output += ('\n')
		self.call_flag = '3:readcalling finished'
	#end of call_pg_lengths()




	def call_length(self):
		"""Function for calling polyG length from called lengths in reads."""
		#set initial values
		self.zygosity = 'unknown'
		self.allel1 = None
		self.allel2 = None
		self.length_dist = {}	#dict to length information
		self.length_dist['total'] = 0	#hold total number of sucessfully called lengths
		length_array = []
		self.pg_lengths_stdev = None
		self.mm_read_count = 0

		#check that read calling has been done and that we have enough sucessfully called reads to work with
		if int(self.call_flag[0]) != 3:
			self.output += 'Warning: Cannot do calling, readcalling is not finished.\n'+self.call_flag+'\n';
			self.call_flag = '6:genotype calling aborted';
			return;
		if  self.sucessfully_called_reads < input.required_read_count:
			self.output += 'Length calling aborted: To few sucessfully called reads.\n';
			self.call_flag = '5:genotype calling aborted';
			return;
		
		#Start genotype calling
		self.call_flag = '4:genotype calling started' #set flag
		#Get distrobution of called lengths, build array for stdev calc and count number of sucessfully called reads with mm in pg
		for read in self.reads:
			if int(read.call_flag[0]) == 0:
				length_array.append(read.called_length)
				self.length_dist['total'] += 1
				try:
					self.length_dist[read.called_length]['count'] += 1
				except KeyError:
					self.length_dist[read.called_length] = {}
					self.length_dist[read.called_length]['count'] = 1
				if read.MM_in_pg: self.mm_read_count += 1

		#calculate the standard deviation of the called lengths population
		from numpy import std
		self.pg_lengths_stdev = std(length_array)#(sum([(float(length) - (sum(length_array) / len(length_array)) ) ** 2 for length in length_array]) / len(length_array)) ** 0.5
		self.output += "Standard deviation of called lengths population:\n"
		self.output += str(self.pg_lengths_stdev)+'\n\n'
		
		#calculate lengthdist width
		self.length_dist_width = max(length_array)-min(length_array)

		#calculate the percentage of successfuly called reads that have a missmatch in the polyG
		self.output += "% successfuly called reads with MM in pG:\n"
		self.percentage_MM = (self.mm_read_count/float(self.length_dist['total']))
		self.output += str(100*self.percentage_MM)+'\n\n'
		
		#calculate percentage support for each allel in distrobution of called lengths
		for length in self.length_dist:
			if length == 'total': continue
			self.length_dist[length]['percentage_support'] = self.length_dist[length]['count']/float(self.length_dist['total'])
		del self.length_dist['total']

		self.allel1_ps = None
		self.allel2_ps = None
		#call homozygous vs heterozygous
		for length in self.length_dist: #check the support for each allel and call if greater than input.percentage_support
			if self.length_dist[length]['percentage_support'] >= input.percentage_support:
				if self.allel1 == None:
					self.allel1 = length
					self.allel1_ps = self.length_dist[length]['percentage_support']
				elif self.allel2 == None:
					self.allel2 = length
					self.allel2_ps = self.length_dist[length]['percentage_support']
				else: raise ValueError, 'cannot call three allels'
		if self.allel1 != None and self.allel2 != None:
			self.zygosity = 'hetrozygous'
		else:
			self.zygosity = 'homozygous'		
			
		#"print" to output
		self.output += "Distribution of called lengths:\n"
		self.output += "len\tcount\t%\n"
		for length in self.length_dist:
			if length == 'total': continue
			self.output += str(length)+'\t'+str(self.length_dist[length]['count'])+'\t'+str(self.length_dist[length]['percentage_support']*100)+'\n'
		self.output += "\n"
		self.output += "Width of called lengths dist:\n"
		self.output += str(self.length_dist_width)+'\n\n'
		self.output += self.pgid+ " is "+ self.zygosity+' allel1='+str(self.allel1)+', allel2='+str(self.allel2)+'\n'
		self.call_flag = '0:genotype calling finished'
	#end of call_length()




	def cleanup(self,level):
		"""Function that removes information that is no longer usefull to save memory.
		levels:
			1. removes bamfile info
			2. removes all read info
			3. removes alignment string only"""
		if level == 1:
			#remove the info from bamfile
			for read in self.reads:
				del read.qname
				del read.flag
				del read.rname
				del read.pos
				del read.mapq
				del read.cigar
				del read.mrnm
				del read.mpos
				del read.isize
				del read.seq
				del read.qual
				del read.opt
				del read.length
		if level == 2:
			#remove all read info
			del self.reads
		if level == 3:
			#remove alignment strings
			for read in self.reads:
				del read.alignment_str
	#end of cleanup()




	def get_per_base_coverage_count(self):
		"""Gets the coverage for each position in polyG"""
		# are done by subprocessing following command:
		# time samtools mpileup data/reference/10-753_WholeGenome/10-753_lanes_12_45_6_S_merged.bam -r 1:17597314-17597321 -f /bubo/nobackup/uppnex/reference/Homo_sapiens/GRCh37/program_files/samtools/chr1.fasta -q 30 # note that its a lot faster if no reference is used!!!
		# make an array or dist with coverage eg. [13, 13, 13, 13, 14, 13] or (4596:13, 4597:13, 4598:13, 4599:13, 4600:14, 4601:13)
		# a dist where we count number of pos with certain coverage eg. (13:5, 14:1) is constructed
		self.per_base_cov = {}
		
		#run subprocess
		samtools = subprocess.Popen(['samtools', 'mpileup', '-r',self.chrom+':'+str(self.pos)+'-'+str( self.pos + self.reference.length -1), input.bamfile.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		pileup_data, errdata = samtools.communicate()
		if samtools.returncode != 0:
			print 'Errorwhile running comand: '+str(' '.join(['samtools', 'mpileup', '-r',self.chrom+':'+str(self.pos)+'-'+str( self.pos + self.reference.length -1), input.bamfile.name])) + ' for polyG='+str(self.pgid)
			print 'Samtools mpileup Error code', samtools.returncode, errdata
			sys.exit()
		pileup_data = StringIO(pileup_data)

		# save coverage for each pos to dist
		line_count = 0
		for line in pileup_data:
			line_count +=1
			[chrom,pos,ref,cov,seq,qual] = line.split('\t')
			cov = int(cov)
			try:
				self.per_base_cov[self.reference.length][cov] += 1
			except KeyError:
				try:
					self.per_base_cov[self.reference.length][cov] = 1
				except KeyError:
					self.per_base_cov[self.reference.length] = {}
					self.per_base_cov[self.reference.length][cov] = 1

		#count positions with zero coverage
		if line_count < self.reference.length:
			try:
				self.per_base_cov[self.reference.length][0] += self.reference.length - line_count
			except KeyError:
				try: self.per_base_cov[self.reference.length][0] = self.reference.length - line_count
				except KeyError:
					self.per_base_cov[self.reference.length] = {}
					self.per_base_cov[self.reference.length][0] = self.reference.length - line_count

		elif line_count > self.reference.length: raise ValueError, 'line count in pilup cant be greater than reflength'
	#end of get_per_base_coverage_count




	def print_output(self):
		"""Function that prints output for the polyG object"""
		#check if first or addition to file
		if type(input.outfile) == str:
			input.outfile2 = input.outfile + '.out.tsv';
			input.outfile = input.outfile + '.perPolygInfo.txt';
			input.outfile = open(input.outfile,'w');
			input.outfile2 = open(input.outfile2,'w');
			input.outfile2.write('# pgid\tstdev_of_called_lengths\tnr_successfully_called_reads\tnr_of_different_called_lengths\twidth_of_called_lengths\tpercentage_MM\tallel1\t%support\tallel2\t%support\n')
		elif type(input.outfile) == file: input.outfile = open(input.outfile.name,'a');  input.outfile2 = open(input.outfile2.name,'a');
		else: raise ValueError, 'input.outfile has to be file or string'

		#write info
		input.outfile.write('###################################################################################################\n')
		if not input.no_bold_bases: input.outfile.write(b_start)
		input.outfile.write(self.pgid+', '+self.chrom+':'+str(self.pos)+', '+str(self.reference.length)+'bp, ' )
		if self.reads_loaded:
			input.outfile.write( str(self.read_count)+' reads loaded' )
			if self.reads_filtered: input.outfile.write( ', '+str(self.filtered_read_count)+' after filtering' )
			if self.call_flag[0] in ['0','3','4','5']: input.outfile.write(', '+str(self.sucessfully_called_reads)+' sucessfully called')
			if not input.no_bold_bases: input.outfile.write(b_end)
			input.outfile.write('.\n')
			if input.print_samtools_command: input.outfile.write( self.get_read_cmd+'\n' )
		if not input.no_bold_bases and not self.reads_loaded: input.outfile.write(b_end);
		input.outfile.write(self.output)
		input.outfile.write('\n')

		input.outfile2.write(self.pgid+'\t'+str(self.pg_lengths_stdev)+'\t'+str(self.sucessfully_called_reads)+'\t'+str(len(self.length_dist))+'\t'+str(self.length_dist_width)+'\t'+str(self.percentage_MM)+'\t'+str(self.allel1)+'\t'+str(self.allel1_ps)+'\t'+str(self.allel2)+'\t'+str(self.allel2_ps)+'\n')

		#close files
		input.outfile.close()
		input.outfile2.close()
	#end of print_output()
#end of PolyG class




class PolyGSummary():
	"""Object holding a summary of called polyG region information"""

	def __init__(self):
		self.summary = ""
		self.per_base_cov = {}
		self.full_pg_read_count = 0
		self.total_pg_lengths = 0
	#end of __init__()




	def add_pre_filter_read_count(self,polyG):
		"""add pre filter read count and total length of polyg regions to summary to be able to calculate coverage during later steps"""
		self.full_pg_read_count += polyG.pre_filter_read_count*polyG.reference.length
		self.total_pg_lengths += polyG.reference.length
	#end of add_pre_filter_read_count()




	def add_per_base_cov(self,polyG):
		"""adds the per base coverage to the summary"""
		for [ref_length, dist] in polyG.per_base_cov.iteritems():
			if ref_length not in self.per_base_cov: self.per_base_cov[ref_length] = {}
			for [cov, count] in dist.iteritems():
				try:
					self.per_base_cov[ref_length][cov] += count
				except KeyError:
					try: self.per_base_cov[ref_length][cov] = count
					except KeyError:
						self.per_base_cov[ref_length] = {}
						self.per_base_cov[ref_length][cov] = count
	#end of add_per_base_cov()




	def summarize(self):
		""""Summarizes all info added to the summary and creates a human readable output (self.summary)"""
		#Per Base and average Coverage
		self.average_coverage = 0
		self.summary = ''
		for ref_length in self.per_base_cov:
			self.summary += 'Per Base Coverage Distribution for pgs with length '+str(ref_length)+':\ncovrage\tpositions\n'
			tempo = self.per_base_cov[ref_length].keys()
			tempo.sort()
			for coverage in tempo:
				count = self.per_base_cov[ref_length][coverage]
				self.average_coverage += coverage*count
				self.summary += str(coverage)+'\t'+str(count)+'\n'
			self.summary += '\n'
		self.average_coverage = float(self.average_coverage)/self.total_pg_lengths
		self.summary += "Average coverage:\n"
		self.summary += str(round(self.average_coverage,2))+'x\n'
		self.summary += '\n'
		
		#Coverage of reads fully covering polyG
		self.summary += "Average coverage of reads covering pg region:\n"
		self.summary += str(round(float(self.full_pg_read_count)/self.total_pg_lengths,2))+'x\n'
	#end of summarize()
#End of summary class




def fasta_get_total_records():
		"""Uses grep and wc to count the total number of fasta records in a fasta file"""
		#grep ">" data/polyGs_found_in_hg19/allInOne2roundNoStars.fa| wc -l , 
		grep = subprocess.Popen(['grep',"pg", input.pgfile.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		total_records, errdata = grep.communicate()
		if grep.returncode != 0:
			print 'Grep Error code', grep.returncode, errdata
			sys.exit()
		wc = subprocess.Popen(['wc',"-l"],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		total_records, errdata = wc.communicate(total_records)
		if wc.returncode != 0:
			print 'Wc Error code', wc.returncode, errdata
			sys.exit()
		return int(total_records.rstrip())




def polyG_generator(pgfile):
	'''Generator for loading info from the pg region infile one at a time.'''
	loaded = 0
	for record in SeqIO.parse(pgfile, "fasta"):
			loaded += 1
			polyG = seqRecord2PolyG(record)
			if input.load_max != 0 and loaded > input.load_max: break
			yield polyG




def seqRecord2PolyG(record):
	"""Function that converts a Bio.SeqRecord object to a polyG object."""
	# Format of fasta entrys
	# >pg9 pg=9|pos=50|len=105|pglen=5|taxid=10090|mol=genomic|chrom=chr1|starposinchrom=13943
	pgid	= record.id
	chrom	= record.description.split('|')[6].split('=')[1][3:]
	pos		= int(record.description.split('|')[7].split('=')[1])
	length	= int(record.description.split('|')[3].split('=')[1])
	seq = str(record.seq)
	return PolyG(pgid,chrom,pos,length,seq)




def bits(i,n): 
	return tuple((0,1)[i>>j & 1] for j in xrange(n-1,-1,-1))

def explainflag(flag):
	""" Function that takes a flag and prints the explanation"""
	########################-Format Description-#########################
	#	Flag	Chr	Description											#
	#1	0x0001	p	the read is paired in sequencing					#
	#2	0x0002	P	the read is mapped in a proper pair					#
	#3	0x0004	u	the query sequence itself is unmapped				#
	#4	0x0008	U	the mate is unmapped								#
	#5	0x0010	r	strand of the query (1 for reverse)					#
	#6	0x0020	R	strand of the mate									#
	#7	0x0040	1	the read is the first read in a pair				#
	#8	0x0080	2	the read is the second read in a pair				#
	#9	0x0100	s	the alignment is not primary						#
	#10	0x0200	f	the read fails platform/vendor quality checks		#
	#11	0x0400	d	the read is either a PCR or an optical duplicate	#
	#####################################################################

	binary = bits(flag,11)
	output = ''
	try:
		if int(binary[-1]):	readtype = "PE"
		else:				readtype = 'SE'
		if int(binary[-2]):	properpair = True
		else:				properpair = False
		if int(binary[-3]):	mapped = False
		else:				mapped = True
		if int(binary[-4]):	matemapped = False
		else:				matemapped = True
		if int(binary[-5]):	strand = 'rev'
		else:				strand = 'fwd'
		if int(binary[-6]):	mate_strand = "rev"
		else:				mate_strand = "fwd"
		if int(binary[-7]):	read1 = True
		else:				read1 = False
		if int(binary[-8]):	read2 = True
		else:				read2 = False
		if int(binary[-9]):	output += 'not primary alignment\t'
		else:				pass
		if int(binary[-10]):output += 'QC failed\t'
		else:				pass
		if int(binary[-11]):output += 'is PCR duplicate'
		else:				pass
	except ValueError:
		output += '.'
	return [readtype, properpair, mapped, matemapped, strand, mate_strand, read1, read2]




#####
#check if run or imported // call main() or not
#####
if __name__ == "__main__":
    main()
#END of script