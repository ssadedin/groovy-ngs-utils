package gngs
import htsjdk.samtools.*
import htsjdk.samtools.util.SequenceUtil
import org.junit.Before
import org.junit.Test

/**
 * Unit tests for DuplexDetector
 */
class DuplexDetectorTest {
    
    DuplexDetector detector
    SAMFileHeader header
    
    @Before
    void setup() {
        detector = new DuplexDetector()
        
        // Create a minimal SAM header
        header = new SAMFileHeader()
        SAMSequenceRecord seq = new SAMSequenceRecord("chr1", 1000000)
        header.addSequence(seq)
    }
    
    /**
     * Helper to create a SAMRecord with specific CIGAR and sequence
     */
    SAMRecord createRecord(String cigar, String bases, String quals = null) {
        SAMRecord record = new SAMRecord(header)
        record.setReadName("test_read")
        record.setReferenceName("chr1")
        record.setAlignmentStart(1000)
        record.setCigarString(cigar)
        record.setReadBases(bases.getBytes())
        
        if (quals == null) {
            // Create default quality scores
            byte[] qualBytes = new byte[bases.length()]
            Arrays.fill(qualBytes, (byte) 30)
            record.setBaseQualities(qualBytes)
        } else {
            record.setBaseQualityString(quals)
        }
        
        record.setMappingQuality(60)
        record.setReadPairedFlag(false)
        record.setReadUnmappedFlag(false)
        
        return record
    }
    
    @Test
    void testExtractSoftClipInfo_noSoftClips() {
        SAMRecord record = createRecord("100M", "A" * 100)
        
        def info = detector.extractSoftClipInfo(record)
        
        assert info.leadingSoftClip == 0
        assert info.trailingSoftClip == 0
        assert info.alignedBases == 100
        assert !info.hasSoftClips()
    }
    
    @Test
    void testExtractSoftClipInfo_trailingSoftClip() {
        SAMRecord record = createRecord("50M50S", "A" * 100)
        
        def info = detector.extractSoftClipInfo(record)
        
        assert info.leadingSoftClip == 0
        assert info.trailingSoftClip == 50
        assert info.alignedBases == 50
        assert info.hasSoftClips()
        assert info.getTotalSoftClip() == 50
    }
    
    @Test
    void testExtractSoftClipInfo_leadingSoftClip() {
        SAMRecord record = createRecord("50S50M", "A" * 100)
        
        def info = detector.extractSoftClipInfo(record)
        
        assert info.leadingSoftClip == 50
        assert info.trailingSoftClip == 0
        assert info.alignedBases == 50
        assert info.hasSoftClips()
    }
    
    @Test
    void testExtractSoftClipInfo_bothEndsSoftClipped() {
        SAMRecord record = createRecord("10S80M10S", "A" * 100)
        
        def info = detector.extractSoftClipInfo(record)
        
        assert info.leadingSoftClip == 10
        assert info.trailingSoftClip == 10
        assert info.alignedBases == 80
        assert info.getTotalSoftClip() == 20
    }
    
    @Test
    void testLengthsMatch_exactMatch() {
        assert detector.lengthsMatch(50, 50)
    }
    
    @Test
    void testLengthsMatch_withinTolerance() {
        assert detector.lengthsMatch(50, 55)  // 10% difference
        assert detector.lengthsMatch(55, 50)
        assert detector.lengthsMatch(100, 95)
        assert detector.lengthsMatch(50, 60)  // 20% difference - now within tolerance
        assert detector.lengthsMatch(60, 50)
    }
    
    @Test
    void testLengthsMatch_outsideTolerance() {
        assert !detector.lengthsMatch(50, 65)  // 30% difference
        assert !detector.lengthsMatch(50, 35)  // 30% difference
    }
    
    @Test
    void testLengthsMatch_zeroValues() {
        assert !detector.lengthsMatch(0, 50)
        assert !detector.lengthsMatch(50, 0)
        assert !detector.lengthsMatch(0, 0)
    }
    
    @Test
    void testExtractSequence() {
        SAMRecord record = createRecord("100M", "ACGTACGTACGT")
        
        assert detector.extractSequence(record, 0, 4) == "ACGT"
        assert detector.extractSequence(record, 4, 8) == "ACGT"
        assert detector.extractSequence(record, 0, 12) == "ACGTACGTACGT"
    }
    
    @Test
    void testExtractSequence_invalidRanges() {
        SAMRecord record = createRecord("100M", "ACGTACGT")
        
        assert detector.extractSequence(record, -1, 4) == ""
        assert detector.extractSequence(record, 0, 100) == ""
        assert detector.extractSequence(record, 5, 4) == ""
    }
    
    @Test
    void testIsDuplexRead_normalReadWithNoSoftClips() {
        SAMRecord record = createRecord("100M", "A" * 100)
        
        assert !detector.isDuplexRead(record)
    }
    
    @Test
    void testIsDuplexRead_softClipTooSmall() {
        // 10bp soft clip, 90bp aligned - not within 20% tolerance
        SAMRecord record = createRecord("90M10S", "A" * 100)
        
        assert !detector.isDuplexRead(record)
    }
    
    @Test
    void testIsDuplexRead_unmappedRead() {
        SAMRecord record = createRecord("*", "A" * 100)
        record.setReadUnmappedFlag(true)
        
        assert !detector.isDuplexRead(record)
    }
    
    @Test
    void testIsDuplexRead_trueDuplexWithMatchingSequences() {
        // Create a sequence where the soft clip is reverse complement of aligned part
        String alignedPart = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"  // 50bp
        String softClipPart = SequenceUtil.reverseComplement(alignedPart)  // 50bp reverse complement
        String fullSeq = alignedPart + softClipPart
        
        SAMRecord record = createRecord("50M50S", fullSeq)
        
        assert detector.isDuplexRead(record)
    }
    
    @Test
    void testIsDuplexRead_sequencesDontMatch() {
        // Soft clip is NOT reverse complement
        String alignedPart = "A" * 50
        String softClipPart = "G" * 50  // Not reverse complement
        String fullSeq = alignedPart + softClipPart
        
        SAMRecord record = createRecord("50M50S", fullSeq)
        
        assert !detector.isDuplexRead(record)
    }
    
    @Test
    void testIsDuplexRead_edgeCaseExactly20PercentDifference() {
        // 50bp aligned, 60bp soft clip (20% difference)
        String alignedPart = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"  // 50bp
        String softClipPart = SequenceUtil.reverseComplement(alignedPart) + "ACGTACGTAC"  // 60bp
        String fullSeq = alignedPart + softClipPart
        
        SAMRecord record = createRecord("50M60S", fullSeq)
        
        boolean result = detector.isDuplexRead(record)
        
        assert result  // Should be detected as duplex at exactly 20% tolerance
    }
    
    @Test
    void testIsDuplexRead_hardClipsAreIgnored() {
        // Hard clips (H) should be ignored, only soft clips (S) matter
        SAMRecord record = createRecord("10H50M50S10H", "A" * 100)
        
        def info = detector.extractSoftClipInfo(record)
        
        assert info.trailingSoftClip == 50
        assert info.alignedBases == 50
        // Hard clips are not counted in soft clip totals
    }
    
    @Test
    void testCustomThresholds() {
        detector.lengthTolerance = 0.20  // 20% tolerance
        detector.alignmentThreshold = 0.70  // 70% match required
        
        assert detector.lengthsMatch(50, 60)  // Now within 20% tolerance
    }
    
    @Test
    void testActualDupexRead() {
        SAM bam = new SAM('src/test/data/example_duplex_reads.small.bam')
        SAMRecord read =  bam.withIterator { i -> i.find { it.readName == '62ec80a8-d4f4-4985-8b9f-e8ec3cae9831'  } }
        assert read != null
        def result = detector.isDuplexRead(read)
        assert result == true : "Read duplex read not detected"
    }
    
    @Test
    void testActualNonDuplexReadNoSoftClips() {
        SAM bam = new SAM('src/test/data/example_duplex_reads.small.bam')
        SAMRecord read =  bam.withIterator { i -> i.find { it.readName == '27a95914-6d76-495b-8618-9c4592140fcd'  } }
        assert read != null

        def result = detector.isDuplexRead(read)
        assert result == false : "Non-duplex read detected as duplex!"
    }

    @Test
    void testActualDuplexRead2() {
        String testReadName = '8e430fd0-7b01-43ef-9f5a-4cc0d2671bf2'  
        SAM bam = new SAM('src/test/data/example_duplex_reads.small.bam')
        SAMRecord read =  bam.withIterator { i -> i.find { it.readName == testReadName} }
        assert read != null

        def result = detector.isDuplexRead(read)
        assert result == true : "Duplex read $testReadName was detected as NOT duplex!"
    }
}
