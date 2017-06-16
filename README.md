# vcf2cytosure

This tool converts a VCF with structural variations to the “.CGH” format used by the
commercial
[CytoSure Interpret Software](https://www.ogt.com/products/246_cytosure_interpret_software)
by OGT (Oxford Gene Technology). CytoSure is made for displaying oligo array measurements.
It works on a set of probes, which this tool emulates. 

Structural variant (SV) types DEL, DUP, TDUP, IDUP, INV, and INS are supported.
For each SVs, at least three probes are generated. If the SV is large enough,
then more probes are generated that are spaced 100 kbp apart.

The SVs are displayed in the following way.

* *DEL* (Deletions): Probes at height -1
* *DUP* (Duplication): Probes at height +1
* *TDUP* (Tandem duplication): Probes at height +1.5
* *IDUP* (Interspersed duplication): Probes at height +0.5
* *INV* (Inversion): Probes at height -0.5
* *INS* (Insertion): Shown as an upwards-pointing "triangle" of 15 probes
* *BND* (Break end): Downwards-pointing triangle

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
Probes are drawn at height 0 if the coverage corresponds to the average
coverage, at +2 for double average coverage, and at -2 for no coverage.
Coverages higher than four times the average coverage are clamped to
a height of +4 and are not shown at their actual height.

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
