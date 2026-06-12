package gngs.tools

import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.atomic.AtomicLong

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMFileWriter
import htsjdk.samtools.SAMRecord

/**
 * Tool to detect and handle Oxford Nanopore duplex reads that were accidentally
 * read in duplex form but processed as simplex reads.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
@Log
class DuplexTrimmer extends ToolBase {
    
    DuplexDetector detector = new DuplexDetector()
    
    // Statistics
    AtomicLong totalReads = new AtomicLong(0)
    AtomicLong duplexReads = new AtomicLong(0)
    AtomicLong rejectedReads = new AtomicLong(0)
    AtomicLong trimmedReads = new AtomicLong(0)
    AtomicLong candidateReads = new AtomicLong(0)
    
    // Trace logging
    String traceReadName = null
    
    // Progress tracking
    static final int PROGRESS_INTERVAL = 5000
    String currentChr = null
    int currentPos = 0
    LinkedBlockingQueue<ReadWrapper> inputQueue = null
    LinkedBlockingQueue<ReadWrapper> outputQueue = null
    
    @Override
    void run() {
        
        int queueSize = opts['queueSize'] ? opts['queueSize'].toInteger() : 10000
        int threads = opts['threads'] ? opts['threads'].toInteger() : Runtime.runtime.availableProcessors()
        
        traceReadName = opts['trace']
        
        // Apply configuration options to detector
        if(opts['lengthTolerance']) {
            detector.lengthTolerance = opts['lengthTolerance'].toDouble()
        }
        
        if(opts['alignmentThreshold']) {
            detector.alignmentThreshold = opts['alignmentThreshold'].toDouble()
        }
        
        log.info "Starting DuplexTrimmer with $threads threads and queue size $queueSize"
        log.info "Length tolerance: ${detector.lengthTolerance}, Alignment threshold: ${detector.alignmentThreshold}"
        
        if(traceReadName) {
            log.info "Trace logging enabled for read: $traceReadName"
        }
        
        this.inputQueue = new LinkedBlockingQueue<>(queueSize)
        this.outputQueue = new LinkedBlockingQueue<>(queueSize)
        
        // Poison pill to signal end
        ReadWrapper POISON = new ReadWrapper(sequenceNumber: -1)
        
        SAM bam = new SAM(opts['i'])
        
        // Parse region if provided
        Region region = null
        if(opts['region']) {
            region = new Region(opts['region'])
            log.info "Processing region: ${region}"
        }
        
        // Start output thread
        Thread outputThread = startOutputThread(outputQueue, POISON, opts['o'], bam)
        
        // Start worker threads
        List<Thread> workers = (1..threads).collect { 
            startWorkerThread(inputQueue, POISON) 
        }
        
        // Input thread (main thread)
        if(region) {
            log.info "Reading input BAM file: ${opts['i']} region: ${region}"
        } else {
            log.info "Reading input BAM file: ${opts['i']}"
        }
        
        long seqNum = 0
        
        Closure processRecord = { SAMRecord record ->
            // Track current position for progress logging (minimal overhead)
            currentChr = record.referenceName
            currentPos = record.alignmentStart
            
            ReadWrapper wrapper = new ReadWrapper(
                record: record,
                sequenceNumber: seqNum++
            )
            
            // Add to both queues immediately - let workers do ALL filtering
            inputQueue.put(wrapper)
            outputQueue.put(wrapper)
            
            long count = totalReads.incrementAndGet()
            
            // Print progress every PROGRESS_INTERVAL reads
            if(count % PROGRESS_INTERVAL == 0) {
                printProgress(count)
            }
        }
        
        if(region) {
            bam.withIterator(region, { i -> i.each(processRecord)})
        } else {
            bam.eachRecord(processRecord)
        }
        
        log.info "Finished reading input. Processed ${totalReads.get()} reads"
        
        // Signal end to workers
        threads.times { inputQueue.put(POISON) }
        
        // Wait for workers to finish
        log.info "Waiting for worker threads to complete..."
        workers*.join()
        
        // Signal end to output thread
        outputQueue.put(POISON)
        outputThread.join()
        
        // Print statistics
        printStatistics()
    }
    
    /**
     * Fast check if record has soft clips without full CIGAR parsing
     */
    @CompileStatic
    boolean hasSoftClips(SAMRecord record) {
        return record.cigar?.cigarElements?.any { 
            it.operator == CigarOperator.SOFT_CLIP 
        }
    }
    
    /**
     * Quick length check without alignment
     */
    @CompileStatic
    boolean lengthsRoughlyMatch(SAMRecord record) {
        def info = detector.extractSoftClipInfo(record)
        if (!info.hasSoftClips()) {
            return false
        }
        
        return detector.lengthsMatch(info.totalSoftClip, info.alignedBases)
    }
    
    /**
     * Start a worker thread that processes reads from the input queue
     */
    Thread startWorkerThread(LinkedBlockingQueue<ReadWrapper> inputQueue, ReadWrapper POISON) {
        Thread.start {
            DuplexDetector threadDetector = new DuplexDetector()
            
            // Copy configuration from main detector
            threadDetector.lengthTolerance = detector.lengthTolerance
            threadDetector.alignmentThreshold = detector.alignmentThreshold
            
            while (true) {
                ReadWrapper wrapper = inputQueue.take()
                
                if (wrapper.sequenceNumber == -1) {
                    inputQueue.put(wrapper)  // Put poison back for other workers
                    break
                }
                
                boolean isTraceRead = traceReadName && wrapper.record.readName == traceReadName
                
                // Fast filter - skip obvious non-candidates
                if (!hasSoftClips(wrapper.record)) {
                    if(isTraceRead) {
                        log.info "TRACE [$traceReadName]: REJECTED - No soft clips found"
                    }
                    wrapper.processed = true
                    continue
                }
                
                if(!lengthsRoughlyMatch(wrapper.record)) {
                    if(isTraceRead) {
                        def info = threadDetector.extractSoftClipInfo(wrapper.record)
                        log.info "TRACE [$traceReadName]: REJECTED - Lengths don't match (soft clip: ${info.totalSoftClip}, aligned: ${info.alignedBases})"
                    }
                    wrapper.processed = true
                    continue
                }
                
                // This is a candidate
                candidateReads.incrementAndGet()
                if(isTraceRead) {
                    log.info "TRACE [$traceReadName]: CANDIDATE - Passed fast filter, performing alignment check"
                }
                
                // Do expensive detection
                wrapper.isDuplex = threadDetector.isDuplexRead(wrapper.record)
                
                if(isTraceRead) {
                    // Calculate alignment score for trace output
                    def info = threadDetector.extractSoftClipInfo(wrapper.record)
                    boolean useTrailing = info.trailingSoftClip >= info.leadingSoftClip
                    int softClipLength = useTrailing ? info.trailingSoftClip : info.leadingSoftClip
                    
                    byte[] bases = wrapper.record.readBases
                    int readLength = bases.length
                    
                    String alignedSeq
                    String softClipSeq
                    
                    if (useTrailing) {
                        int alignedEnd = readLength - softClipLength
                        int alignedStart = Math.max(0, alignedEnd - 50)
                        alignedSeq = threadDetector.extractSequence(wrapper.record, alignedStart, alignedEnd)
                        
                        int softClipStart = alignedEnd
                        int softClipEnd = Math.min(readLength, softClipStart + 50)
                        softClipSeq = threadDetector.extractSequence(wrapper.record, softClipStart, softClipEnd)
                    } else {
                        int softClipEnd = softClipLength
                        int softClipStart = Math.max(0, softClipEnd - 50)
                        softClipSeq = threadDetector.extractSequence(wrapper.record, softClipStart, softClipEnd)
                        
                        int alignedStart = softClipEnd
                        int alignedEnd = Math.min(readLength, alignedStart + 50)
                        alignedSeq = threadDetector.extractSequence(wrapper.record, alignedStart, alignedEnd)
                    }
                    
                    // Reverse complement the soft clip sequence for visual comparison
                    String softClipRC = htsjdk.samtools.util.SequenceUtil.reverseComplement(softClipSeq)
                    
                    double alignmentScore = threadDetector.alignSequences(alignedSeq, softClipSeq)
                    
                    log.info "TRACE [$traceReadName]: Alignment comparison (aligned vs RC of soft clip):"
                    log.info "TRACE [$traceReadName]: Aligned:  $alignedSeq"
                    log.info "TRACE [$traceReadName]: SoftClip: $softClipRC"
                    
                    if(wrapper.isDuplex) {
                        log.info "TRACE [$traceReadName]: DUPLEX DETECTED - Alignment score: ${String.format('%.3f', alignmentScore)} (threshold: ${threadDetector.alignmentThreshold})"
                    } else {
                        log.info "TRACE [$traceReadName]: NOT DUPLEX - Alignment score: ${String.format('%.3f', alignmentScore)} (threshold: ${threadDetector.alignmentThreshold})"
                    }
                }
                
                if (wrapper.isDuplex) {
                    duplexReads.incrementAndGet()
                    // Determine which end is the duplex side
                    def duplexInfo = threadDetector.extractSoftClipInfo(wrapper.record)
                    wrapper.duplexOnTrailingEnd = (duplexInfo.trailingSoftClip >= duplexInfo.leadingSoftClip)
                    wrapper.record = applyAction(wrapper.record, wrapper.duplexOnTrailingEnd)
                    
                    if(isTraceRead) {
                        String action = opts['action'] ?: 'reject'
                        log.info "TRACE [$traceReadName]: Action '$action' applied"
                    }
                }
                
                wrapper.processed = true
            }
        }
    }
    
    /**
     * Start the output thread that writes processed reads in order from the output queue
     */
    Thread startOutputThread(LinkedBlockingQueue<ReadWrapper> outputQueue, 
                            ReadWrapper POISON, String outputFile, SAM bam) {
        Thread.start {
            log.info "Writing output to: $outputFile"
            
            bam.withWriter(outputFile, true) { SAMFileWriter writer ->
                
                while (true) {
                    ReadWrapper wrapper = outputQueue.take()
                    
                    if (wrapper.sequenceNumber == -1) {
                        break
                    }
                    
                    // Wait until this specific record is processed
                    while (!wrapper.processed) {
                        Thread.sleep(1)
                    }
                    
                    // Write immediately in queue order - queue order IS the correct order
                    // null record means it was rejected
                    if (wrapper.record != null) {
                        writer.addAlignment(wrapper.record)
                    }
                }
            }
            
            log.info "Finished writing output"
        }
    }
    
    /**
     * Apply the configured action to a duplex read
     */
    @CompileStatic
    SAMRecord applyAction(SAMRecord record, boolean duplexOnTrailingEnd) {
        String action = opts['action'] ?: 'reject'
        
        switch(action) {
            case 'reject':
                rejectedReads.incrementAndGet()
                return null  // Signal to skip this read
                
            case 'trim':
                trimmedReads.incrementAndGet()
                return trimDuplexEnd(record, duplexOnTrailingEnd)
                
            case 'keep':
                return record  // Keep as-is
                
            default:
                log.warning "Unknown action: $action, defaulting to reject"
                rejectedReads.incrementAndGet()
                return null
        }
    }
    
    /**
     * Remove the soft clip from only the end identified as the duplex portion
     */
    @CompileStatic
    SAMRecord trimDuplexEnd(SAMRecord record, boolean duplexOnTrailingEnd) {
        def info = detector.extractSoftClipInfo(record)
        
        List cigarElements = record.cigar.cigarElements
        
        // Only remove the soft clip on the duplex end
        List newCigar
        int trimStart = 0
        int trimEnd = record.readLength
        
        if (duplexOnTrailingEnd) {
            // Remove trailing soft clip only
            newCigar = cigarElements.findAll { el ->
                // Keep all elements except the last one if it's a soft clip
                el != cigarElements.last() || el.operator != CigarOperator.SOFT_CLIP
            }
            trimEnd = record.readLength - info.trailingSoftClip
        } else {
            // Remove leading soft clip only
            newCigar = cigarElements.findAll { el ->
                // Keep all elements except the first one if it's a soft clip
                el != cigarElements.first() || el.operator != CigarOperator.SOFT_CLIP
            }
            trimStart = info.leadingSoftClip
        }
        
        if (newCigar.isEmpty()) {
            return null  // Nothing left after trimming
        }
        
        // Trim the bases and qualities
        if (trimStart < trimEnd) {
            byte[] newBases = Arrays.copyOfRange(record.readBases, trimStart, trimEnd)
            byte[] newQuals = Arrays.copyOfRange(record.baseQualities, trimStart, trimEnd)
            
            record.readBases = newBases
            record.baseQualities = newQuals
        }
        
        // Update the CIGAR after trimming bases
        record.cigar = new htsjdk.samtools.Cigar(newCigar)
        
        return record
    }
    
    /**
     * Print progress statistics
     */
    void printProgress(long count) {
        long candidates = candidateReads.get()
        long duplex = duplexReads.get()
        long rejected = rejectedReads.get()
        long trimmed = trimmedReads.get()
        
        double candidatePercent = count > 0 ? 100.0 * candidates / count : 0.0
        double duplexPercent = candidates > 0 ? 100.0 * duplex / candidates : 0.0
        
        String position = currentChr ? "${currentChr}:${currentPos}" : "N/A"
        
        int inputQueueSize = inputQueue?.size() ?: 0
        int outputQueueSize = outputQueue?.size() ?: 0
        int queueCapacity = inputQueue?.remainingCapacity() ?: 0
        int totalCapacity = inputQueueSize + queueCapacity
        
        double inputUtilization = totalCapacity > 0 ? 100.0 * inputQueueSize / totalCapacity : 0.0
        double outputUtilization = totalCapacity > 0 ? 100.0 * outputQueueSize / totalCapacity : 0.0
        
        log.info String.format(
            "Progress: %,d reads at %s | Queue: in=%d/%d (%.1f%%) out=%d/%d (%.1f%%) | %,d candidates (%.2f%%) | %,d duplex (%.2f%%) | %,d rejected | %,d trimmed",
            count, position, inputQueueSize, totalCapacity, inputUtilization, outputQueueSize, totalCapacity, outputUtilization,
            candidates, candidatePercent, duplex, duplexPercent, rejected, trimmed
        )
    }
    
    /**
     * Print statistics about the run
     */
    void printStatistics() {
        System.err.println ""
        System.err.println "=" * 80
        System.err.println "Duplex Trimmer Statistics"
        System.err.println "=" * 80
        System.err.println "Total reads processed:    ${totalReads.get()}"
        System.err.println "Candidate reads:          ${candidateReads.get()} (${String.format('%.2f', 100.0 * candidateReads.get() / totalReads.get())}%)"
        System.err.println "Duplex reads detected:    ${duplexReads.get()} (${String.format('%.2f', 100.0 * duplexReads.get() / totalReads.get())}%)"
        System.err.println "Reads rejected:           ${rejectedReads.get()}"
        System.err.println "Reads trimmed:            ${trimmedReads.get()}"
        System.err.println "=" * 80
        System.err.println ""
    }
    
    static void main(String[] args) {
        
        Utils.setupAcceleratedDeflaters()
        
        cli('Detect and handle Oxford Nanopore duplex reads', args) {
            i 'Input BAM/SAM file', args: 1, required: true
            o 'Output BAM file', args: 1, required: true
            action 'Action to take for duplex reads: reject, trim, or keep', args: 1, required: false
            threads 'Number of worker threads', args: 1, required: false
            queueSize 'Size of processing queue', args: 1, required: false
            lengthTolerance 'Length tolerance for duplex detection (default 0.10)', args: 1, required: false
            alignmentThreshold 'Alignment threshold for duplex detection (default 0.40)', args: 1, required: false
            region 'Region to process (format: chr:start-end)', args: 1, required: false
            trace 'Read name to trace through processing', args: 1, required: false
        }
    }
}

/**
 * Wrapper class to hold a SAMRecord and its processing state
 */
@CompileStatic
class ReadWrapper {
    SAMRecord record
    long sequenceNumber
    volatile boolean processed = false
    volatile boolean isDuplex = false
    volatile boolean skipProcessing = false
    volatile boolean duplexOnTrailingEnd = false
}
