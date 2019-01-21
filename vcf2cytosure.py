#!/usr/bin/env python3
"""
Convert structural variants in a VCF to CGH (CytoSure) format
"""

import argparse
import sys
import logging
import gzip
import math
from collections import namedtuple, defaultdict
from io import StringIO
from lxml import etree
from cyvcf2 import VCF


from constants import *

__version__ = '0.4.3'

logger = logging.getLogger(__name__)

Event = namedtuple('Event', ['chrom', 'start', 'end', 'type', 'info'])

def events(variants):
	"""Iterate over variants and yield Events"""

	for variant in variants:
		if len(variant.ALT) != 1:
			continue
		chrom = variant.CHROM
		if chrom not in CONTIG_LENGTHS:
			continue
		start = variant.start
		sv_type = variant.INFO.get('SVTYPE')
		if variant.INFO.get("END"):
			end = variant.INFO.get('END')
			if start >= end:
				tmp=int(end)
				end=start
				start=tmp

			logger.debug('%s at %s:%s-%s (%s bp)', sv_type, chrom, start+1, end, end - start)
			assert len(variant.REF) == 1

			yield Event(chrom=chrom, start=start, end=end, type=sv_type, info=dict(variant.INFO))
		else:

			if ":" in variant.ALT[0] and ("[" in  variant.ALT[0] or "]" in  variant.ALT[0]):
				chrom2=variant.ALT[0].split(":")[0].split("[")[-1].split("]")[-1]
			else:
		
				print ("invalid variant  type {}: skipping".format(variant.ALT[0]))

			if chrom2 != chrom:
				logger.debug('%s at %s:%s', sv_type, chrom, start+1)
				yield Event( chrom=chrom, start=start, end=None, type=sv_type, info=dict(variant.INFO) )

			else:
				end=int(variant.ALT[0].split(":")[1].split("[")[0].split("]")[0])
				if start >= end:
					tmp=int(end)
					end=start
					start=tmp

				logger.debug('%s at %s:%s', sv_type, chrom, start+1)
				yield Event( chrom=chrom, start=start, end=end, type=sv_type, info=dict(variant.INFO) )

def strip_template(path):
	"""
	Read in the template CGH file and strip it of everything that we don't need.

	Return the lxml.etree object.
	"""
	tree = etree.parse(path)

	# Remove all aberrations
	parent = tree.xpath('/data/cgh/submission')[0]
	for aberration in parent.xpath('aberration'):
		parent.remove(aberration)

	# Remove all except the first probe (in the order in which they occur in
	# the file) on each chromosome. Chromosomes without probes are not
	# clickable in the CytoSure UI.
	parent = tree.xpath('/data/cgh/probes')[0]
	seen = set()
	for probe in parent:
		chrom = probe.attrib.get('chromosome')
		if not chrom or chrom in seen:
			parent.remove(probe)
		else:
			seen.add(chrom)

	# Remove all segments
	parent = tree.xpath('/data/cgh/segmentation')[0]
	for segment in parent:
		parent.remove(segment)

	return tree


def make_probe(parent, chromosome, start, end, height, text):
	probe = etree.SubElement(parent, 'probe')
	probe.attrib.update({
		'name': text,
		'chromosome': CHROM_RENAME.get(chromosome, chromosome),
		'start': str(start + 1),
		'stop': str(end),
		'normalized': '{:.3f}'.format(-height),
		'smoothed': '0.0',
		'smoothed_normalized': '-0.25',
		'sequence': 'AACCGGTT',
	})

	red = 1000
	green = red * 2**height

	spot = etree.SubElement(probe, 'spot')
	spot.attrib.update({
		'index': '1',
		'row': '1',
		'column': '1',
		'red': str(red),
		'green': '{:.3f}'.format(green),
		'gSNR': '100.0',
		'rSNR': '100.0',
		'outlier': 'false',
	})
	return probe


def make_segment(parent, chromosome, start, end, height):
	segment = etree.SubElement(parent, 'segment')
	segment.attrib.update({
		'chrId': CHROM_RENAME.get(chromosome, chromosome),
		'numProbes': '100',
		'start': str(start + 1),
		'stop': str(end),
		'average': '{:.3f}'.format(-height),  # CytoSure inverts the sign
	})
	return segment


def make_aberration(parent, chromosome, start, end, comment=None, method='converted from VCF',
		confirmation=None, n_probes=0, copy_number=99):
	"""
	comment -- string
	method -- short string
	confirmation -- string
	"""
	aberration = etree.SubElement(parent, 'aberration')
	aberration.attrib.update(dict(
		chr=CHROM_RENAME.get(chromosome, chromosome),
		start=str(start + 1),
		stop=str(end),
		maxStart=str(start + 1),
		maxStop=str(end),
		copyNumber=str(copy_number),
		initialClassification='Unclassified',
		finalClassification='Unclassified',
		inheritance='Not_tested',
		numProbes=str(n_probes),
		startProbe='',
		stopProbe='',
		maxStartProbe='',
		maxStopProbe='',

		# TODO fill in the following values with something sensible
		automationLevel='1.0',
		baseline='0.0',
		mosaicism='0.0',
		gain='true',
		inheritanceCoverage='0.0',
		logRatio='-0.4444',  # mean log ratio
		method=method,
		p='0.003333',  # p-value
		sd='0.2222',  # standard deviation
	))
	if comment:
		e = etree.SubElement(aberration, 'comments')
		e.text = comment
	if confirmation:
		e = etree.SubElement(aberration, 'confirmation')
		e.text = confirmation
	return aberration


def spaced_probes(start, end, probe_spacing=PROBE_SPACING):
	"""
	Yield nicely spaced positions along the interval (start, end).
	- start and end are always included
	- at least three positions are included
	"""
	l = end - start
	n = l // probe_spacing
	spacing = l / max(n, 2)  # float division
	i = 0
	pos = start
	while pos <= end:
		yield pos
		i += 1
		pos = start + int(i * spacing)


def probe_point(center, height=2.5, width=5001, steps=15):
	"""
	Yield (pos, height) pairs that "draw" a triangular shape (pointing upwards)
	"""
	pos_step = (width - 1) // (steps - 1)
	height_step = height / ((steps - 1) // 2)
	for i in range(-(steps // 2), steps // 2 + 1):
		yield center + i * pos_step, height - height_step * abs(i) + 0.1


def format_comment(info):
	comment = ''
	for k, v in sorted(info.items()):
		if k in ('CSQ', 'SVTYPE'):
			continue
		comment += '\n{}: {}'.format(k, v)
	return comment


def merge_intervals(intervals):
	"""Merge overlapping intervals into a single one"""
	events = [(coord[0], 'START') for coord in intervals]
	events.extend((coord[1], 'STOP') for coord in intervals)
	events.sort()
	active = 0
	start = 0
	for pos, what in events:
		# Note adjacent 'touching' events are merged because 'START' < 'STOP'
		if what == 'START':
			if active == 0:
				start = pos
			active += 1
		else:
			active -= 1
			if active == 0:
				yield (start, pos)


def complement_intervals(intervals, chromosome_length):
	"""
	>>> list(complement_intervals([(0, 1), (3, 4), (18, 20)], 20))
	[(1, 3), (4, 18)]
	"""
	prev_end = 0
	for start, end in intervals:
		if prev_end != start:
			yield prev_end, start
		prev_end = end
	if prev_end != chromosome_length:
		yield prev_end, chromosome_length


def add_probes_between_events(probes, chr_intervals):
	for chrom, intervals in chr_intervals.items():
		if chrom not in CONTIG_LENGTHS:
			continue
		intervals = merge_intervals(intervals)
		for start, end in complement_intervals(intervals, CONTIG_LENGTHS[chrom]):
			for pos in spaced_probes(start, end, probe_spacing=200000):
				# CytoSure does not display probes at height=0.0
				make_probe(probes, chrom, pos, pos + 60, 0.01, 'between events')


class CoverageRecord:
	__slots__ = ('chrom', 'start', 'end', 'coverage')

	def __init__(self, chrom, start, end, coverage):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.coverage = coverage


def parse_coverages(path):
	with open(path) as f:
		for line in f:
			if line.startswith('#'):
				continue
			content=line.split('\t')
			chrom=content[0]
			start=content[1]
			end=content[2]
			coverage=content[3]
			start = int(start)
			end = int(end)
			coverage = float(coverage)
			yield CoverageRecord(chrom, start, end, coverage)

def retrieve_snp(content,args):
	snp_data=[]
	snp_data.append(content[0])
	snp_data.append(int(content[1]))
	snp_data.append(float( content[7].split(";{}=".format(args.dp))[-1].split(";")[0] ))

	return (snp_data)

def parse_snv_coverages(args):
		snv_list=[]
		if args.snv.endswith(".gz"):
			for line  in gzip.open(args.snv):
				if line[0] == "#" or not ";{}=".format(args.dp) in line:
					continue
				content=line.strip().split()
				snv_list.append(retrieve_snp(content,args))

		elif args.snv.endswith(".vcf"):
			for line in open(args.snv):
				if line[0] == "#" or not ";{}=".format(args.dp) in line:
					continue
				content=line.strip().split()
				snv_list.append(retrieve_snp(content,args))

		else:
			print ("only .vcf or gziped vcf is allowed, exiting")

		for snv in snv_list:
			chrom = snv[0]
			start = snv[1]
			end = start+1
			coverage = snv[2]
			yield CoverageRecord(chrom, start, end, coverage)

def group_by_chromosome(records):
	"""
	Group records by their .chrom attribute.

	Yield pairs (chromosome, list_of_records) where list_of_records
	are the consecutive records sharing the same chromosome.
	"""
	prev_chrom = None
	chromosome_records = []
	for record in records:
		if record.chrom != prev_chrom:
			if chromosome_records:
				yield prev_chrom, chromosome_records
				chromosome_records = []
		chromosome_records.append(record)
		prev_chrom = record.chrom
	if chromosome_records:
		yield prev_chrom, chromosome_records


def bin_coverages(coverages, n):
	"""
	Reduce the number of coverage records by re-binning
	each *n* coverage values into a new single bin.

	The coverages are assumed to be from a single chromosome.
	"""
	chrom = coverages[0].chrom
	for i in range(0, len(coverages), n):
		records = coverages[i:i+n]
		cov = sum(r.coverage for r in records) / len(records)
		yield CoverageRecord(chrom,	records[0].start, records[-1].end, cov)


def subtract_intervals(records, intervals):
	"""
	Yield only those records that fall outside of the given intervals.
	"""
	events = [(r.start, 'rec', r) for r in records]
	events.extend((i[0], 'istart', None) for i in intervals)
	events.extend((i[1], 'iend', None) for i in intervals)
	events.sort()
	inside = False
	for pos, typ, record in events:
		if typ == 'istart':
			inside = True
		elif typ == 'iend':
			inside = False
		elif not inside:
			yield record


def add_coverage_probes(probes, path,args):
	"""
	probes -- <probes> element
	path -- path to tab-separated file with coverages
	"""
	logger.info('Reading %r ...', path)
	if args.coverage:
		coverages = [r for r in parse_coverages(path) if r.chrom in CONTIG_LENGTHS]
	else:
		coverages = [r for r in parse_snv_coverages(args) if r.chrom in CONTIG_LENGTHS]

	mean_coverage = sum(r.coverage for r in coverages) / len(coverages)
	logger.info('Mean coverage is %.2f', mean_coverage)

	n = 0
	for chromosome, records in group_by_chromosome(coverages):
		coverage_factor = 2
		if args.sex == 'male' and ( chromosome == 'Y' or chromosome == 'X'):
			coverage_factor = 1

		n_intervals = N_INTERVALS[chromosome]
		for record in subtract_intervals(bin_coverages(records,args.bins), n_intervals):
			height = min(coverage_factor * record.coverage / mean_coverage - 2, MAX_HEIGHT)
			if height == 0.0:
				height = 0.01
			make_probe(probes, record.chrom, record.start, record.end, height, 'coverage')
			n += 1
	logger.info('Added %s coverage probes', n)


#apply filtering
def variant_filter(variants, min_size=5000,max_frequency=0.01, frequency_tag='FRQ'):

	for variant in variants:

		end = variant.INFO.get('END')
		if end and not variant.INFO.get('SVTYPE') == 'TRA':

			if abs( int(end) - variant.start) <= min_size:
				# Too short
				continue

		elif variant.INFO.get('SVTYPE') == 'BND':
			bnd_chrom, bnd_pos = variant.ALT[0][2:-1].split(':')

			bnd_pos = int(variant.ALT[0].split(':')[1].split("]")[0].split("[")[0])
			bnd_chrom= variant.ALT[0].split(':')[0].split("]")[-1].split("[")[-1]

			if bnd_chrom == variant.CHROM and abs(bnd_pos - variant.start) < min_size:
				continue

		elif variant.INFO.get('SVTYPE') == 'TRA':

			bnd_pos = variant.INFO.get('END')
			bnd_chrom =variant.INFO.get('CHR2');

		frequency = variant.INFO.get(frequency_tag)
		if frequency is not None and frequency > max_frequency:
			continue

		yield variant

#retrieve the sample id, assuming single sample vcf
def retrieve_sample_id(vcf_path):
	sample=vcf_path.split("/")[-1].split("_")[0].split(".")[0]
	return(sample)


def main():
	logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	parser = argparse.ArgumentParser("VCF2cytosure - convert SV vcf files to cytosure")

	group = parser.add_argument_group('Filtering')
	group.add_argument('--size', default=1000, type=int,help='Minimum variant size. Default: %(default)s')
	group.add_argument('--frequency', default=0.01, type=float,help='Maximum frequency. Default: %(default)s')
	group.add_argument('--frequency_tag', default='FRQ', type=str,help='Frequency tag of the info field. Default: %(default)s')
	group.add_argument('--no-filter', dest='do_filtering', action='store_false',default=True,help='Disable any filtering')

	group = parser.add_argument_group('Input')
	group.add_argument('--coverage',help='Coverage file')
	group.add_argument('--sex',required=False, default='female', help='Sample sex male/female. Default: %(default)s')
	group.add_argument('--vcf',required=True,help='VCF file')
	group.add_argument('--bins',type=int,default=20,help='the number of coverage bins per probes default=20')
	group.add_argument('--snv',type=str,help='snv vcf file, use coverage annotation to position the height of the probes(cannot be used together with --coverage)')
	group.add_argument('--dp',type=str,default="DP",help='read depth tag of snv vcf file,this option is only  used if you use snv to set the heigth of the probes. The dp tag is a tag which is used to retrieve the depth of coverage across the snv (default=DP)')
	group.add_argument('--maxbnd',type=int,default=10000,help='Maxixmum BND size, BND events exceeding this size are discarded')
	group.add_argument('--out',help='output file (default = the prefix of the input vcf)')

	group.add_argument('-V','--version',action='version',version="%(prog)s "+__version__ ,
			   help='Print program version and exit.')
	# parser.add_argument('xml', help='CytoSure design file')
	args= parser.parse_args()

	if args.coverage and args.snv:
		print ("--coverage and --snv cannot be combined")
		quit()

	if not args.out:
		args.out=".".join(args.vcf.split(".")[0:len(args.vcf.split("."))-1])+".cgh"
	parser = etree.XMLParser(remove_blank_text=True)

	sex_male = "false"
	promega_sex = 'Female'
	if args.sex == "male":
		sex_male = 'true'
		promega_sex = 'Male'

	sample_id=retrieve_sample_id(args.vcf)
	tree = etree.parse(StringIO(CGH_TEMPLATE.format(sample_id,sex_male,promega_sex,sex_male)), parser)

	segmentation = tree.xpath('/data/cgh/segmentation')[0]
	probes = tree.xpath('/data/cgh/probes')[0]
	submission = tree.xpath('/data/cgh/submission')[0]

	chr_intervals = defaultdict(list)
	vcf = VCF(args.vcf)
	if args.do_filtering:
		vcf = variant_filter(vcf,min_size=args.size,max_frequency=args.frequency,frequency_tag=args.frequency_tag)
	n = 0
	for event in events(vcf):
		height = ABERRATION_HEIGHTS[event.type]
		end = event.end
		make_segment(segmentation, event.chrom, event.start, end, height)
		comment = format_comment(event.info)
		if "rankScore" in event.info:
			rank_score = int(event.info['RankScore'].partition(':')[2])
		else:
			rank_score =0

		#occ=0
		#if args.frequency_tag in event.info:
		#	occ=event.info[args.frequency_tag]
		occ=0
		if "OCC" in event.info:
			occ=event.info["OCC"]

		if event.type in ("INV",'INS', 'BND',"TRA") and not event.end:
			continue
			#pass
		elif event.type in ("INV",'INS', 'BND',"TRA") and (abs(event.start-event.end) > args.maxbnd ):
			#pass
			continue

		make_aberration(submission, event.chrom, event.start, end, confirmation=event.type,
			comment=comment, n_probes=occ, copy_number=rank_score)


		
		chr_intervals[event.chrom].append((event.start, event.end))
		# show probes at slightly different height than segments
		for pos in spaced_probes(event.start, event.end - 1):
			make_probe(probes, event.chrom, pos, pos + 60, height, event.type)
		n += 1
	if args.coverage or args.snv:
		add_coverage_probes(probes, args.coverage, args)
	else:
		add_probes_between_events(probes, chr_intervals)

	tree.write(args.out, pretty_print=True)
	logger.info('Wrote %d variants to CGH', n)


if __name__ == '__main__':
	main()
