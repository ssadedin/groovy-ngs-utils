
def cols = ['Func','Gene','ExonicFunc','AAChange','Conserved','SegDup','ESP5400_ALL','g1000','dbSNP132','AVSIFT','LJB_PhyloP','LJB_PhyloP_Pred','LJB_SIFT','LJB_SIFT_Pred','LJB_PolyPhen2','LJB_PolyPhen2_Pred','LJB_LRT','LJB_LRT_Pred','LJB_MutationTaster','LJB_MutationTaster_Pred','LJB_GERP++','Chr','Start','End','Ref','Obs','Otherinfo']

def annovar = new CSV("153566_S2_L001_R1_001.fastq.trim.atrim.reorder.realign.recal.target.merge.vcf.av.refgene.exome_summary.csv",
    columnNames: cols,
    readFirstLine: false)


Iterator i = annovar.iterator()

line1 = i.next()
println "Line 1 is $line1"

line2 = i.next()
println "Line 2 is $line2"
