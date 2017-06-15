#!/usr/bin/env python3
"""
Convert structural variants in a VCF to CGH (CytoSure) format
"""

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import sys
import logging
import math
from collections import namedtuple, defaultdict
from io import StringIO
from lxml import etree
from cyvcf2 import VCF

PROBE_SPACING = 100000
MAX_HEIGHT = 4

logger = logging.getLogger(__name__)


# GRCh37
CONTIG_LENGTHS = {
	'1':  249_250_621,
	'2':  243_199_373,
	'3':  198_022_430,
	'4':  191_154_276,
	'5':  180_915_260,
	'6':  171_115_067,
	'7':  159_138_663,
	'8':  146_364_022,
	'9':  141_213_431,
	'10': 135_534_747,
	'11': 135_006_516,
	'12': 133_851_895,
	'13': 115_169_878,
	'14': 107_349_540,
	'15': 102_531_392,
	'16':  90_354_753,
	'17':  81_195_210,
	'18':  78_077_248,
	'19':  59_128_983,
	'20':  63_025_520,
	'21':  48_129_895,
	'22':  51_304_566,
	'X':  155_270_560,
	'Y':   59_373_566,
}

Event = namedtuple('Event', ['chrom', 'start', 'end', 'type', 'info'])

CGH_TEMPLATE = """
<data formatVersion="2">
<pgdData><pgdDataEntry key="SPECIMEN_TYPE" value="BLASTOMERE"/></pgdData>
<noResults>
</noResults>
<pgd_reagents/>
<cgh mother="-1" father="-1" genomeBuild="hg19" softwareVersion="4.8.32" batched="false">
  <submission design="031035" feFile="dummy.txt" cghFile="dummy.cgh" scanDate="1462520414000" barcode="253103511677_1_3" sampleCy3="true">
  <notes/>
  <sample sampleId="2016-08276" male="false"><phenotype/></sample>
  <reference sampleId="Promega Female" male="false"><phenotype/></reference>
  <extra>
    <datum category="Nanodrop" type="Sample DNA (ng)" dataType="Float"/>
    <datum category="Sample Extraction" type="Sample Arrival Date" dataType="Date"/>
    <datum category="Sample Extraction" type="Specimen Type" dataType="List">Blood</datum>
    <datum category="Labelling" type="Lab User" dataType="String"/>
    <datum category="Hyb &amp; Wash" type="Cot1 Batch" dataType="String"/>
    <datum category="Sample Extraction" type="Extraction Method" dataType="String"/>
    <datum category="General" type="Reference Concentration" dataType="Float"/>
    <datum category="Sample Extraction" type="A260/A280" dataType="Float"/>
    <datum category="General" type="Assigned Technologist" dataType="String">MF</datum>
    <datum category="Nanodrop" type="Reference DNA (pmoles)" dataType="Float"/>
    <datum category="Labelling" type="Columns Used" dataType="List">Qiagen</datum>
    <datum category="Nanodrop" type="Sample DNA (A260/A280)" dataType="Float"/>
    <datum category="Labelling" type="Column Batch No" dataType="String"/>
    <datum category="Sample Extraction" type="Extracted By" dataType="String"/>
    <datum category="Hyb &amp; Wash" type="Hyb Protocol" dataType="String"/>
    <datum category="Labelling" type="Experiment Date" dataType="Date"/>
    <datum category="General" type="Case Status" dataType="List"/>
    <datum category="Labelling" type="Lab Notes" dataType="Text"/>
    <datum category="Hyb &amp; Wash" type="Hyb Buffer Batch" dataType="String"/>
    <datum category="Nanodrop" type="Reference DNA (A260/A280)" dataType="Float"/>
    <datum category="Labelling" type="Labelling Reagent Batch No" dataType="String"/>
    <datum category="Hyb &amp; Wash" type="Wash Protocol" dataType="String"/>
    <datum category="Labelling" type="Labelling Protocol" dataType="List">Enzo</datum>
    <datum category="Nanodrop" type="Reference DNA (ng)" dataType="Float"/>
    <datum category="Nanodrop" type="Sample DNA (pmoles)" dataType="Float"/>
  </extra>
</submission>
<excludedRegions>
</excludedRegions>
<qc>
<aqf key="SPIKES" value="null"/>
<aqf key="DLRSPREAD" value="0.1"/>
<aqf key="RED_SIGNAL_INTENSITY" value="1000.0"/>
<aqf key="GREEN_SIGNAL_INTENSITY" value="1000.0"/>
<aqf key="BG_NOISE_RED" value="5.0"/>
<aqf key="BG_NOISE_GREEN" value="5.0"/>
<aqf key="RED_SNR" value="100.0"/>
<aqf key="GREEN_SNR" value="100.0"/>
<aqf key="COMPARATIVE_SIGNAL_INTENSITY" value="0.11111111111"/>
<aqf key="REPRODUCIBILITY_GREEN" value="0.1111111"/>
<aqf key="REPRODUCIBILITY_RED" value="0.1111111"/>
<aqf key="NEG_CONTROL_RED" value="1.1111111"/>
<aqf key="NEG_CONTROL_GREEN" value="1.1111111"/>
<aqf key="FLAG_PERC" value="0.00722363"/>
<aqf key="SAT_PERC" value="0.0"/>
<aqf key="WAVINESS" value="0.0011111"/>
<aqf key="SNP_TROUGH_PEAK_RATIO" value="null"/>
<aqf key="SNP_GREY_AREA" value="null"/>
<aqf key="SNP_Red_Signal_Intensity" value="NaN"/>
<aqf key="SNP_Green_Signal_Intensity" value="NaN"/>
<aqf key="SNP_Ratio_Separation" value="0.0"/>
<aqf key="Percentage_Homozygosity" value="NaN"/>
<aqf key="SD" value="0.111111"/>
</qc>
<probes>
</probes>
<segmentation type="NORMALIZED"></segmentation>
</cgh>
</data>
"""


class HelpfulArgumentParser(ArgumentParser):
	"""An ArgumentParser that prints full help on errors."""

	def __init__(self, *args, **kwargs):
		if 'formatter_class' not in kwargs:
			kwargs['formatter_class'] = RawDescriptionHelpFormatter
		super().__init__(*args, **kwargs)

	def error(self, message):
		self.print_help(sys.stderr)
		args = {'prog': self.prog, 'message': message}
		self.exit(2, '%(prog)s: error: %(message)s\n' % args)


def parse_vcf(path):
	for variant in VCF(path):
		if len(variant.ALT) != 1:
			continue
		chrom = variant.CHROM
		start = variant.start
		sv_type = variant.INFO.get('SVTYPE')
		if sv_type in ('DEL', 'DUP', 'IDUP', 'TDUP', 'INV', 'INS'):
			end = variant.INFO.get('END')
			assert start <= end
			if sv_type not in ('IDUP', 'INS'):
				assert sv_type == variant.ALT[0][1:-1]  # <DUP> etc.
			logger.info('%s at %s:%s-%s (%s bp)', sv_type, chrom, start+1, end, end - start)
			assert len(variant.REF) == 1

			yield Event(chrom=chrom, start=start, end=end, type=sv_type, info=dict(variant.INFO))
		#elif vtype == 'BND':  # Breakend
			#print(variant.orientation, variant.chr, variant.connectingSequence, variant.pos, variant.remoteOrientation, variant.withinMainAssembly)
			#assert r.INFO.get('SVTYPE') == 'BND'
			#break
		#else:
			##print(dir(variant))
			#keys = frozenset(['END', 'IMPRECISE', 'SVTYPE', 'SVLEN'])
			#print({k:variant.INFO.get(k) for k in keys})
			#if variant.INFO.get('SVTYPE') != vtype:
				#print(variant, vtype)
				#assert False


def strip_template(path):
	"""
	Read in the template CGH file and strip it of everything that we don’t need.

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
		'chromosome': chromosome,
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
		'chrId': chromosome,
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
		chr=chromosome,
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


def triangle_probes(center, height=2.5, width=5001, steps=15):
	"""
	Yield (pos, height) pairs that "draw" a triangular shape (pointing upwards)
	"""
	pos_step = (width - 1) // (steps - 1)
	height_step = height / ((steps - 1) // 2)
	for i in range(-(steps // 2), steps // 2 + 1):
		yield center + i * pos_step, height - height_step * abs(i) + 0.1


def format_comment(info: dict) -> str:
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


def main():
	logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	parser = HelpfulArgumentParser(description=__doc__)
	parser.add_argument('--coverage',
		help='Coverage file')
	parser.add_argument('vcf',
		help="VCF file")
	# parser.add_argument('xml', help='CytoSure design file')
	args = parser.parse_args()

	parser = etree.XMLParser(remove_blank_text=True)
	tree = etree.parse(StringIO(CGH_TEMPLATE), parser)
	segmentation = tree.xpath('/data/cgh/segmentation')[0]
	probes = tree.xpath('/data/cgh/probes')[0]
	submission = tree.xpath('/data/cgh/submission')[0]

	chr_intervals = defaultdict(list)
	for event in parse_vcf(args.vcf):
		height = {'DEL': -1.0, 'DUP': +1.0, 'TDUP': +1.5, 'IDUP': +0.5, 'INV': -0.5, 'INS': -1.5}[event.type]
		make_segment(segmentation, event.chrom, event.start, event.end, height)

		comment = format_comment(event.info)
		rank_score = int(event.info['RankScore'].partition(':')[2])
		make_aberration(submission, event.chrom, event.start, event.end, confirmation=event.type,
			comment=comment, n_probes=event.info['OCC'], copy_number=rank_score)

		if event.type == 'INS':
			for pos, height in triangle_probes(event.start):
				make_probe(probes, event.chrom, pos, pos + 60, height, event.type)
		else:
			chr_intervals[event.chrom].append((event.start, event.end))
			# show probes at slightly different height than segments
			height *= 1.05
			for pos in spaced_probes(event.start, event.end - 1):
				make_probe(probes, event.chrom, pos, pos + 60, height, event.type)

	if args.coverage:
		i = 0
		for line in open(args.coverage):
			if line.startswith('#'):
				continue
			i += 1
			if i % 10 != 0:
				continue
			chrom, start, end, coverage, _ = line.split('\t')
			start = int(start) + 1
			end = int(end)
			center = (end + start) // 2
			coverage = float(coverage)
			height = min(2 * coverage / 30 - 2, MAX_HEIGHT)  # FIXME hardcoded avg cov
			make_probe(probes, chrom, center - 30, center + 30, height, 'coverage')

	else:
		def probes_between_events(chrom, start, end):
			for pos in spaced_probes(start, end, probe_spacing=200000):
				# CytoSure does not display probes at height=0.0
				make_probe(probes, chrom, pos, pos + 60, 0.01, 'between events')

		for chrom, intervals in chr_intervals.items():
			if chrom not in CONTIG_LENGTHS:
				continue
			intervals = merge_intervals(intervals)
			prev_end = 1
			for start, end in intervals:
				probes_between_events(chrom, prev_end, start - 1)
				prev_end = end
			probes_between_events(chrom, prev_end, CONTIG_LENGTHS[chrom])

	tree.write(sys.stdout.buffer, pretty_print=True)


if __name__ == '__main__':
	main()
