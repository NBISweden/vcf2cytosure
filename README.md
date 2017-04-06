# vcf2cytosure

Convert VCF with structural variations to CytoSure format



## Notes

- CGH is XML. The company does not have/does not give out schema.
- The parser in CytoSure is quite robust. It is fine, for example, to re-format
  the XML with

      xmllint --format file.cgh > pretty.xml
- one probe can be made up of more than one spot:

    <probe ...>
        <spot ... />
        <spot ... />
    </probe>


- There are `<probe>` elements without sequence, without chromosome name:

      <probe name="probename" sequence="">
        <spot index="56774" row="334" column="164" red="123.4" green="345.6" gSNR="1.0" rSNR="50.0" outlier="false"/>
      </probe>

- coordinates are 1-based (since stop - start = 59 and length of sequence is 60)
- Removing the 'sequence' attribute does not work, but setting it to a fake one does
- log2 ratio is computed as log2(green/red)
- normalization segments:
      <segment chrId="1" numProbes="10" start="7000" stop="1000000" average="-0.003"/>
- segments can overlap each other
