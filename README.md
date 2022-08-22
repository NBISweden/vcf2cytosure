# vcf2cytosure

This tool converts a VCF with structural variations to the “.CGH” format used by the
commercial
[CytoSure Interpret Software](https://www.ogt.com/products/246_cytosure_interpret_software)
by OGT (Oxford Gene Technology). CytoSure is made for displaying oligo array measurements.
It works on a set of probes, which this tool emulates. 

## Usage
Cytosure requires an input vcf file. Optionally, a coverage bed or snp vcf file could be used to visualise the coverage across the genome.

    INSTALL:
        git  clone https://github.com/NBISweden/vcf2cytosure.git
        cd vcf2cytosure
        pip install -e .
	
    or use Singularity  
    
    	singularity pull shub://J35P312/vcf2cytosure

    RUNNING:
        
    vcf2cytosure --vcf <input.vcf> --out <output.cgh>
    
    optionally

    vcf2cytosure --vcf <input.vcf> --out <output.cgh> --coverge <overage.bed> --sex <male|female>

    or:

    vcf2cytosure --vcf <input.vcf> --out <output.cgh> --snv <snv.vcf> --sex <male|female>

The coverage bed file may be created using TIDDIT(https://github.com/J35P312/TIDDIT)

    TIDDIT --cov -b <input.bam> -o <coverage_pefix>

The binning of the input coverage file may be controlled using the --bins parameter:

	vcf2cytosure --vcf <input.vcf> --out <output.cgh> --snv <snv.vcf> --bins 50

Here 50 coverage bins will be pooled into one probe. The number of probes affect the amount of detail and resolution in the analysis.
A large number of probes will make cytosure sluggish.

## Structural variants
Structural variant (SV) types DEL, DUP, TDUP, IDUP, INV, and INS are supported.
For each SVs, at least three probes are generated. If the SV is large enough,
then more probes are generated that are spaced 100 kbp apart.

The SVs are displayed in the following way.

* *DEL* (Deletions): Probes at height -1
* *DUP* (Duplication): Probes at height +1
* *TDUP* (Tandem duplication): Probes at height +1.5

For each SV, an *aberration* record is also generated. The attributes are
filled in the following way:

* *“Confirmation Method”* is set to the type of the SV (DEL, INS, DUP, etc.)
* *“# Probes”* is set to the OCC value found in the input INFO field (number
  of occurrences in the reference population)
* The *“Comment”* field contains all the INFO attributes that have not
  otherwise been used for one of the above fields.
* The copyNumber field is set to the RankScore value from the VCF INFO field
* The other attributes are set to some arbitrary, but constant value.

If the `--coverage` option is used, probe heights will represent coverage.
Coverage is represented as log2 ratios of input bin coverage relative to 
all bin coverages, limited at [MIN_HEIGHT, MAX_HEIGHT]. 
This means probes are drawn at height 0 if the coverage corresponds to the average
coverage, at 0.58 for a heterozygote duplication, 1 for double average coverage, -1 for a heterozygous deletion
and at around -MAX_HEIGHT for little or no coverage.

High or low average coverages are clamped to a height of +/-4 and are not shown at their actual height.

If the `--coverage` option is not used, evenly spaced probes with height 0.01
(height 0.0 would not be shown by CytoSure) are generated for the areas between
SVs.



## Notes on the file format

- CGH is in XML format. The company does not seem to have a schema file.
- The parser in CytoSure accepts reformatted XML files. This is ok:

      xmllint --format file.cgh > pretty.xml
- One probe can be made up of more than one spot:

      <probe ...>
        <spot ... />
        <spot ... />
      </probe>

- The probes in the file do not have to be sorted by chromosome and/or
  coordinate.
- There are `<probe>` elements without sequence, without chromosome name:

      <probe name="probename" sequence="">
        <spot index="56774" row="334" column="164" red="123.4" green="345.6" gSNR="1.0" rSNR="50.0" outlier="false"/>
      </probe>

- Coordinates are 1-based (since stop - start = 59 and length of sequence is 60)
- Removing the 'sequence' attribute does not work, but setting it to a fake one does
- The log2 ratio is computed as log2(green/red)
- Normalization segments look like this:

      <segment chrId="1" numProbes="10" start="7000" stop="1000000" average="-0.003"/>
- Segments can overlap each other
- Names "X" and "Y" are not used for the corresponding chromosomes. Instead,
  the `X` chromosome is stored as `23`, they `Y` chromosome as `24` (in probes,
  segments and aberrations).

- N_INTERVALS for hg38 were taken from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz.
