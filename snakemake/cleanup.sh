UNPOLISHED=unpolished
mkdir -p $UNPOLISHED

ln -rs unpolished.fasta ${UNPOLISHED}/unpolished.fasta
ln -rs unpolished.fasta.fai ${UNPOLISHED}/unpolished.fasta.fai

mv unpolished.ref_reads.pbmm2.bam ${UNPOLISHED}/unpolished.ref.bam
mv unpolished.ref_reads.pbmm2.bam.bai ${UNPOLISHED}/unpolished.ref.bam.bai

mv unpolished.ref_reads.polish_ref.fasta ${UNPOLISHED}/refpolished.fasta

mv unpolished.ref_reads.polish_ref.query_reads.bwa.bam ${UNPOLISHED}/refpolished.query.bam
mv unpolished.ref_reads.polish_ref.query_reads.bwa.bam.bai ${UNPOLISHED}/refpolished.query.bam.bai

# --------------------------------------------------------

CONSENSUS=consensus
mkdir -p $CONSENSUS

rm consensus.query_reads.longranger.vcf.gz
rm consensus.ref_reads.longshot.vcf.gz.tbi
rm consensus.query_reads.longranger.bam
rm consensus.query_reads.longranger.bam.bai
rm consensus.ref_reads.pbmm2.bam
rm consensus.ref_reads.pbmm2.bam.bai
rm high_confidence_hets.vcf.gz
rm high_confidence_hets.vcf.gz.tbi

mv consensus.fasta ${CONSENSUS}/consensus.fasta
mv consensus.fasta.fai ${CONSENSUS}/consensus.fasta.fai

mv consensus.query_reads.variants.bam ${CONSENSUS}/consensus.query.bam
mv consensus.query_reads.variants.bam.bai ${CONSENSUS}/consensus.query.bam.bai
mv consensus.query_reads.variants.vcf.gz ${CONSENSUS}/consensus.query.vcf.gz

mv consensus.ref_reads.variants.bam ${CONSENSUS}/consensus.ref.bam
mv consensus.ref_reads.variants.bam.bai ${CONSENSUS}/consensus.ref.bam.bai
mv consensus.ref_reads.variants.vcf.gz ${CONSENSUS}/consensus.ref.vcf.gz

mv ref_reads_haplotag.bam ${CONSENSUS}/ref.haplotag.bam
mv ref_reads_haplotag.bam.bai ${CONSENSUS}/ref.haplotag.bam.bai
mv query_reads_haplotag.bam ${CONSENSUS}/query.haplotag.bam
mv query_reads_haplotag.bam.bai ${CONSENSUS}/query.haplotag.bam.bai

mv high_confidence_hets.whatshap.vcf.gz ${CONSENSUS}/high_confidence_hets.phased.vcf.gz
mv high_confidence_hets.whatshap.vcf.gz.tbi ${CONSENSUS}/high_confidence_hets.phased.vcf.gz.tbi

PILONC=${CONSENSUS}/pilon
mkdir -p $PILONC

mv unpolished.ref_reads.polish_ref_query_reads_pilonBadCoverage.wig ${PILONC}/refpolished.BadCoverage.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilon.changes ${PILONC}/refpolished.changes
mv unpolished.ref_reads.polish_ref_query_reads_pilonChanges.wig ${PILONC}/refpolished.Changes.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilonClippedAlignments.wig ${PILONC}/refpolished.ClippedAlignments.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilonCopyNumber.wig ${PILONC}/refpolished.CopyNumber.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilonCoverage.wig ${PILONC}/refpolished.Coverage.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilonDeltaCoverage.wig ${PILONC}/refpolished.DeltaCoverage.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilonDipCoverage.wig ${PILONC}/refpolished.DipCoverage.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilonGC.wig ${PILONC}/refpolished.GC.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilonPctBad.wig ${PILONC}/refpolished.PctBad.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilonPhysicalCoverage.wig ${PILONC}/refpolished.PhysicalCoverage.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilonPilon.bed ${PILONC}/refpolished.Pilon.bed
mv unpolished.ref_reads.polish_ref_query_reads_pilonUnconfirmed.wig ${PILONC}/refpolished.Unconfirmed.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilonWeightedMq.wig ${PILONC}/refpolished.WeightedMq.wig
mv unpolished.ref_reads.polish_ref_query_reads_pilonWeightedQual.wig ${PILONC}/refpolished.WeightedQual.wig


# --------------------------------------------------------

READS=reads
mkdir -p $READS

mv query_reads_haplotype1.bam ${READS}/query_hap1.bam
mv query_reads_haplotype1.bai ${READS}/query_hap1.bam.bai 
mv query_reads_haplotype2.bam ${READS}/query_hap2.bam
mv query_reads_haplotype2.bam.bai ${READS}/query_hap2.bam.bai 
mv query_reads_unassigned.bam ${READS}/query_unassigned.bam
mv query_reads_unassigned.bam.bai ${READS}/query_unassigned.bam.bai

mv ref_reads_haplotype1.bam ${READS}/ref_hap1.bam
mv ref_reads_haplotype1.bam.bai ${READS}/ref_hap1.bam.bai 
mv ref_reads_haplotype2.bam ${READS}/ref_hap2.bam
mv ref_reads_haplotype2.bam.bai ${READS}/ref_hap2.bam.bai 
mv ref_reads_unassigned.bam ${READS}/ref_unassigned.bam 
mv ref_reads_unassigned.bam.bai ${READS}/ref_unassigned.bam.bai

# --------------------------------------------------------

ITER0=iter0
mkdir -p $ITER0

ln -rs ${CONSENSUS}/consensus.fasta ${ITER0}/unpolished.hap1.fasta
ln -rs ${CONSENSUS}/consensus.fasta ${ITER0}/unpolished.hap2.fasta

mv iter0.hap1.unpolished.ref_reads_haplotype1_iter0.pbmm2.bam ${ITER0}/unpolished.ref.hap1.bam
mv iter0.hap1.unpolished.ref_reads_haplotype1_iter0.pbmm2.bam.bai ${ITER0}/unpolished.ref.hap1.bam.bai
mv iter0.hap2.unpolished.ref_reads_haplotype2_iter0.pbmm2.bam ${ITER0}/unpolished.ref.hap2.bam
mv iter0.hap2.unpolished.ref_reads_haplotype2_iter0.pbmm2.bam.bai ${ITER0}/unpolished.ref.hap2.bam.bai

mv iter0.hap1.unpolished.ref_reads_haplotype1_iter0.hap_polish_ref.fasta ${ITER0}/refpolished.hap1.fasta
mv iter0.hap2.unpolished.ref_reads_haplotype2_iter0.hap_polish_ref.fasta ${ITER0}/refpolished.hap2.fasta

mv iter0.hap1.unpolished.ref_reads_haplotype1_iter0.hap_polish_ref.query_reads_haplotype1_iter0.bwa.bam ${ITER0}/refpolished.query.hap1.bam
mv iter0.hap1.unpolished.ref_reads_haplotype1_iter0.hap_polish_ref.query_reads_haplotype1_iter0.bwa.bam.bai ${ITER0}/refpolished.query.hap1.bam.bai
mv iter0.hap2.unpolished.ref_reads_haplotype2_iter0.hap_polish_ref.query_reads_haplotype2_iter0.bwa.bam ${ITER0}/refpolished.query.hap2.bam
mv iter0.hap2.unpolished.ref_reads_haplotype2_iter0.hap_polish_ref.query_reads_haplotype2_iter0.bwa.bam.bai ${ITER0}/refpolished.query.hap2.bam.bai

mv iter0.hap1.consensus.fasta ${ITER0}/consensus.hap1.fasta
mv iter0.hap1.consensus.fasta.fai ${ITER0}/consensus.hap1.fasta.fai
mv iter0.hap2.consensus.fasta ${ITER0}/consensus.hap2.fasta
mv iter0.hap2.consensus.fasta.fai ${ITER0}/consensus.hap2.fasta.fai

PILON0=${ITER0}/pilon
mkdir -p $PILON0

mv *haplotype1_iter0_pilonBadCoverage.wig ${PILON0}/hap1.BadCoverage.wig
mv *haplotype1_iter0_pilon.changes ${PILON0}/hap1.changes
mv *haplotype1_iter0_pilonChanges.wig ${PILON0}/hap1.Changes.wig
mv *haplotype1_iter0_pilonClippedAlignments.wig ${PILON0}/hap1.ClippedAlignments.wig
mv *haplotype1_iter0_pilonCopyNumber.wig ${PILON0}/hap1.CopyNumber.wig
mv *haplotype1_iter0_pilonCoverage.wig ${PILON0}/hap1.Coverage.wig
mv *haplotype1_iter0_pilonDeltaCoverage.wig ${PILON0}/hap1.DeltaCoverage.wig
mv *haplotype1_iter0_pilonDipCoverage.wig ${PILON0}/hap1.DipCoverage.wig
mv *haplotype1_iter0_pilonGC.wig ${PILON0}/hap1.GC.wig
mv *haplotype1_iter0_pilonPctBad.wig ${PILON0}/hap1.PctBad.wig
mv *haplotype1_iter0_pilonPhysicalCoverage.wig ${PILON0}/hap1.PhysicalCoverage.wig
mv *haplotype1_iter0_pilonPilon.bed ${PILON0}/hap1.Pilon.bed
mv *haplotype1_iter0_pilonUnconfirmed.wig ${PILON0}/hap1.Unconfirmed.wig
mv *haplotype1_iter0_pilonWeightedMq.wig ${PILON0}/hap1.WeightedMq.wig
mv *haplotype1_iter0_pilonWeightedQual.wig ${PILON0}/hap1.WeightedQual.wig
mv *haplotype2_iter0_pilonBadCoverage.wig ${PILON0}/hap2.BadCoverage.wig
mv *haplotype2_iter0_pilon.changes ${PILON0}/hap2.changes
mv *haplotype2_iter0_pilonChanges.wig ${PILON0}/hap2.Changes.wig
mv *haplotype2_iter0_pilonClippedAlignments.wig ${PILON0}/hap2.ClippedAlignments.wig
mv *haplotype2_iter0_pilonCopyNumber.wig ${PILON0}/hap2.CopyNumber.wig
mv *haplotype2_iter0_pilonCoverage.wig ${PILON0}/hap2.Coverage.wig
mv *haplotype2_iter0_pilonDeltaCoverage.wig ${PILON0}/hap2.DeltaCoverage.wig
mv *haplotype2_iter0_pilonDipCoverage.wig ${PILON0}/hap2.DipCoverage.wig
mv *haplotype2_iter0_pilonGC.wig ${PILON0}/hap2.GC.wig
mv *haplotype2_iter0_pilonPctBad.wig ${PILON0}/hap2.PctBad.wig
mv *haplotype2_iter0_pilonPhysicalCoverage.wig ${PILON0}/hap2.PhysicalCoverage.wig
mv *haplotype2_iter0_pilonPilon.bed ${PILON0}/hap2.Pilon.bed
mv *haplotype2_iter0_pilonUnconfirmed.wig ${PILON0}/hap2.Unconfirmed.wig
mv *haplotype2_iter0_pilonWeightedMq.wig ${PILON0}/hap2.WeightedMq.wig
mv *haplotype2_iter0_pilonWeightedQual.wig ${PILON0}/hap2.WeightedQual.wig


# --------------------------------------------------------

ITER1=iter1
mkdir -p $ITER1

ln -rs ${ITER1}/consensus.hap1.fasta ${ITER1}/unpolished.hap1.fasta
ln -rs ${ITER1}/consensus.hap2.fasta ${ITER1}/unpolished.hap2.fasta

mv iter1.hap1.unpolished.ref_reads_haplotype1_iter1.pbmm2.bam ${ITER1}/unpolished.ref.hap1.bam
mv iter1.hap1.unpolished.ref_reads_haplotype1_iter1.pbmm2.bam.bai ${ITER1}/unpolished.ref.hap1.bam.bai
mv iter1.hap2.unpolished.ref_reads_haplotype2_iter1.pbmm2.bam ${ITER1}/unpolished.ref.hap2.bam
mv iter1.hap2.unpolished.ref_reads_haplotype2_iter1.pbmm2.bam.bai ${ITER1}/unpolished.ref.hap2.bam.bai

mv iter1.hap1.unpolished.ref_reads_haplotype1_iter1.hap_polish_ref.fasta ${ITER1}/refpolished.hap1.fasta
mv iter1.hap2.unpolished.ref_reads_haplotype2_iter1.hap_polish_ref.fasta ${ITER1}/refpolished.hap2.fasta

mv iter1.hap1.unpolished.ref_reads_haplotype1_iter1.hap_polish_ref.query_reads_haplotype1_iter1.bwa.bam ${ITER1}/refpolished.query.hap1.bam
mv iter1.hap1.unpolished.ref_reads_haplotype1_iter1.hap_polish_ref.query_reads_haplotype1_iter1.bwa.bam.bai ${ITER1}/refpolished.query.hap1.bam.bai
mv iter1.hap2.unpolished.ref_reads_haplotype2_iter1.hap_polish_ref.query_reads_haplotype2_iter1.bwa.bam ${ITER1}/refpolished.query.hap2.bam
mv iter1.hap2.unpolished.ref_reads_haplotype2_iter1.hap_polish_ref.query_reads_haplotype2_iter1.bwa.bam.bai ${ITER1}/refpolished.query.hap2.bam.bai

mv iter1.hap1.consensus.fasta ${ITER1}/consensus.hap1.fasta
mv iter1.hap1.consensus.fasta.fai ${ITER1}/consensus.hap1.fasta.fai
mv iter1.hap2.consensus.fasta ${ITER1}/consensus.hap2.fasta
mv iter1.hap2.consensus.fasta.fai ${ITER1}/consensus.hap2.fasta.fai

PILON1=${ITER1}/pilon
mkdir -p $PILON1

mv *haplotype1_iter1_pilonBadCoverage.wig ${PILON1}/hap1.BadCoverage.wig
mv *haplotype1_iter1_pilon.changes ${PILON1}/hap1.changes
mv *haplotype1_iter1_pilonChanges.wig ${PILON1}/hap1.Changes.wig
mv *haplotype1_iter1_pilonClippedAlignments.wig ${PILON1}/hap1.ClippedAlignments.wig
mv *haplotype1_iter1_pilonCopyNumber.wig ${PILON1}/hap1.CopyNumber.wig
mv *haplotype1_iter1_pilonCoverage.wig ${PILON1}/hap1.Coverage.wig
mv *haplotype1_iter1_pilonDeltaCoverage.wig ${PILON1}/hap1.DeltaCoverage.wig
mv *haplotype1_iter1_pilonDipCoverage.wig ${PILON1}/hap1.DipCoverage.wig
mv *haplotype1_iter1_pilonGC.wig ${PILON1}/hap1.GC.wig
mv *haplotype1_iter1_pilonPctBad.wig ${PILON1}/hap1.PctBad.wig
mv *haplotype1_iter1_pilonPhysicalCoverage.wig ${PILON1}/hap1.PhysicalCoverage.wig
mv *haplotype1_iter1_pilonPilon.bed ${PILON1}/hap1.Pilon.bed
mv *haplotype1_iter1_pilonUnconfirmed.wig ${PILON1}/hap1.Unconfirmed.wig
mv *haplotype1_iter1_pilonWeightedMq.wig ${PILON1}/hap1.WeightedMq.wig
mv *haplotype1_iter1_pilonWeightedQual.wig ${PILON1}/hap1.WeightedQual.wig
mv *haplotype2_iter1_pilonBadCoverage.wig ${PILON1}/hap2.BadCoverage.wig
mv *haplotype2_iter1_pilon.changes ${PILON1}/hap2.changes
mv *haplotype2_iter1_pilonChanges.wig ${PILON1}/hap2.Changes.wig
mv *haplotype2_iter1_pilonClippedAlignments.wig ${PILON1}/hap2.ClippedAlignments.wig
mv *haplotype2_iter1_pilonCopyNumber.wig ${PILON1}/hap2.CopyNumber.wig
mv *haplotype2_iter1_pilonCoverage.wig ${PILON1}/hap2.Coverage.wig
mv *haplotype2_iter1_pilonDeltaCoverage.wig ${PILON1}/hap2.DeltaCoverage.wig
mv *haplotype2_iter1_pilonDipCoverage.wig ${PILON1}/hap2.DipCoverage.wig
mv *haplotype2_iter1_pilonGC.wig ${PILON1}/hap2.GC.wig
mv *haplotype2_iter1_pilonPctBad.wig ${PILON1}/hap2.PctBad.wig
mv *haplotype2_iter1_pilonPhysicalCoverage.wig ${PILON1}/hap2.PhysicalCoverage.wig
mv *haplotype2_iter1_pilonPilon.bed ${PILON1}/hap2.Pilon.bed
mv *haplotype2_iter1_pilonUnconfirmed.wig ${PILON1}/hap2.Unconfirmed.wig
mv *haplotype2_iter1_pilonWeightedMq.wig ${PILON1}/hap2.WeightedMq.wig
mv *haplotype2_iter1_pilonWeightedQual.wig ${PILON1}/hap2.WeightedQual.wig

# --------------------------------------------------------

ITER2=iter2
mkdir -p $ITER2

ln -rs ${ITER2}/consensus.hap1.fasta ${ITER2}/unpolished.hap1.fasta
ln -rs ${ITER2}/consensus.hap2.fasta ${ITER2}/unpolished.hap2.fasta

mv iter2.hap1.unpolished.ref_reads_haplotype1_iter2.pbmm2.bam ${ITER2}/unpolished.ref.hap1.bam
mv iter2.hap1.unpolished.ref_reads_haplotype1_iter2.pbmm2.bam.bai ${ITER2}/unpolished.ref.hap1.bam.bai
mv iter2.hap2.unpolished.ref_reads_haplotype2_iter2.pbmm2.bam ${ITER2}/unpolished.ref.hap2.bam
mv iter2.hap2.unpolished.ref_reads_haplotype2_iter2.pbmm2.bam.bai ${ITER2}/unpolished.ref.hap2.bam.bai

mv iter2.hap1.unpolished.ref_reads_haplotype1_iter2.hap_polish_ref.fasta ${ITER2}/refpolished.hap1.fasta
mv iter2.hap2.unpolished.ref_reads_haplotype2_iter2.hap_polish_ref.fasta ${ITER2}/refpolished.hap2.fasta

mv iter2.hap1.unpolished.ref_reads_haplotype1_iter2.hap_polish_ref.query_reads_haplotype1_iter2.bwa.bam ${ITER2}/refpolished.query.hap1.bam
mv iter2.hap1.unpolished.ref_reads_haplotype1_iter2.hap_polish_ref.query_reads_haplotype1_iter2.bwa.bam.bai ${ITER2}/refpolished.query.hap1.bam.bai
mv iter2.hap2.unpolished.ref_reads_haplotype2_iter2.hap_polish_ref.query_reads_haplotype2_iter2.bwa.bam ${ITER2}/refpolished.query.hap2.bam
mv iter2.hap2.unpolished.ref_reads_haplotype2_iter2.hap_polish_ref.query_reads_haplotype2_iter2.bwa.bam.bai ${ITER2}/refpolished.query.hap2.bam.bai

mv iter2.hap1.consensus.fasta ${ITER2}/consensus.hap1.fasta
mv iter2.hap1.consensus.fasta.fai ${ITER2}/consensus.hap1.fasta.fai
mv iter2.hap2.consensus.fasta ${ITER2}/consensus.hap2.fasta
mv iter2.hap2.consensus.fasta.fai ${ITER2}/consensus.hap2.fasta.fai

PILON2=${ITER2}/pilon
mkdir -p $PILON2

mv *haplotype1_iter2_pilonBadCoverage.wig ${PILON2}/hap1.BadCoverage.wig
mv *haplotype1_iter2_pilon.changes ${PILON2}/hap1.changes
mv *haplotype1_iter2_pilonChanges.wig ${PILON2}/hap1.Changes.wig
mv *haplotype1_iter2_pilonClippedAlignments.wig ${PILON2}/hap1.ClippedAlignments.wig
mv *haplotype1_iter2_pilonCopyNumber.wig ${PILON2}/hap1.CopyNumber.wig
mv *haplotype1_iter2_pilonCoverage.wig ${PILON2}/hap1.Coverage.wig
mv *haplotype1_iter2_pilonDeltaCoverage.wig ${PILON2}/hap1.DeltaCoverage.wig
mv *haplotype1_iter2_pilonDipCoverage.wig ${PILON2}/hap1.DipCoverage.wig
mv *haplotype1_iter2_pilonGC.wig ${PILON2}/hap1.GC.wig
mv *haplotype1_iter2_pilonPctBad.wig ${PILON2}/hap1.PctBad.wig
mv *haplotype1_iter2_pilonPhysicalCoverage.wig ${PILON2}/hap1.PhysicalCoverage.wig
mv *haplotype1_iter2_pilonPilon.bed ${PILON2}/hap1.Pilon.bed
mv *haplotype1_iter2_pilonUnconfirmed.wig ${PILON2}/hap1.Unconfirmed.wig
mv *haplotype1_iter2_pilonWeightedMq.wig ${PILON2}/hap1.WeightedMq.wig
mv *haplotype1_iter2_pilonWeightedQual.wig ${PILON2}/hap1.WeightedQual.wig
mv *haplotype2_iter2_pilonBadCoverage.wig ${PILON2}/hap2.BadCoverage.wig
mv *haplotype2_iter2_pilon.changes ${PILON2}/hap2.changes
mv *haplotype2_iter2_pilonChanges.wig ${PILON2}/hap2.Changes.wig
mv *haplotype2_iter2_pilonClippedAlignments.wig ${PILON2}/hap2.ClippedAlignments.wig
mv *haplotype2_iter2_pilonCopyNumber.wig ${PILON2}/hap2.CopyNumber.wig
mv *haplotype2_iter2_pilonCoverage.wig ${PILON2}/hap2.Coverage.wig
mv *haplotype2_iter2_pilonDeltaCoverage.wig ${PILON2}/hap2.DeltaCoverage.wig
mv *haplotype2_iter2_pilonDipCoverage.wig ${PILON2}/hap2.DipCoverage.wig
mv *haplotype2_iter2_pilonGC.wig ${PILON2}/hap2.GC.wig
mv *haplotype2_iter2_pilonPctBad.wig ${PILON2}/hap2.PctBad.wig
mv *haplotype2_iter2_pilonPhysicalCoverage.wig ${PILON2}/hap2.PhysicalCoverage.wig
mv *haplotype2_iter2_pilonPilon.bed ${PILON2}/hap2.Pilon.bed
mv *haplotype2_iter2_pilonUnconfirmed.wig ${PILON2}/hap2.Unconfirmed.wig
mv *haplotype2_iter2_pilonWeightedMq.wig ${PILON2}/hap2.WeightedMq.wig
mv *haplotype2_iter2_pilonWeightedQual.wig ${PILON2}/hap2.WeightedQual.wig




# --------------------------------------------------------

rm iter*.hap*.unpolished.fasta
rm iter*.hap*.unpolished.fasta.fai
rm query_reads_haplotype*_iter*.fastq.gz

rm query_reads_haplotype*.bam
rm query_reads_haplotype*.bam.bai


# --------------------------------------------------------



