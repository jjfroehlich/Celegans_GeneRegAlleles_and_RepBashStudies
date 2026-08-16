<b>"<i>C. elegans</i> Gene Regulatory Alleles and Reporter Bashing Studies"</b>
<br> Froehlich, JJ; Rajewsky, N (2023). microPublication Biology. <a href="https://doi.org/10.17912/micropub.biology.000709" target="_blank" rel="noopener noreferrer">10.17912/micropub.biology.000709</a>.

<br> [R script](https://github.com/jjfroehlich/Celegans_GeneRegAlleles_and_RepBashStudies/blob/main/Froehlich_Celegans_GeneRegAlleles_and_RepBashStudies.R) to replicate analyses and figures.
<br> [HTML notebook with code and figures](https://jjfroehlich.github.io/Celegans_GeneRegAlleles_and_RepBashStudies/)

<br><br>
<b>Main steps of script:</b>
<ul>
  <li>Downloads “classical alleles” and genome/transcript annotations from WormBase (WS284)</li>
  <li>Cleans transcript models to conservatively remove overlaps (for example UTRs overlapping CDS) and generates clean, non-overlapping genomic features: CDS, UTRs, introns, and intergenic regions</li>
  <li>Classifies "classical alleles" by overlap with these features</li>
  <li>Exports <code>.bed</code> files and an Excel spreadsheet for downstream curation. Then re-loads manually curated infos</li>
  <li>Generates all plots shown in the publication, and additional summary plots as <code>.pdf</code></li>
  <li>Creates automated GViz browser snapshots of alleles with surrounding annotation to verify</li>
</ul>
