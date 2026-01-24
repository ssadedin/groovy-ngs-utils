package gngs.tools

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMFileWriter
import htsjdk.samtools.SAMRecord

import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.atomic.AtomicLong

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
        
        LinkedBlockingQueue<ReadWrapper> inputQueue = new LinkedBlockingQueue<>(queueSize)
        LinkedBlockingQueue<ReadWrapper> outputQueue = new LinkedBlockingQueue<>(queueSize)
        
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
            boolean isTraceRead = traceReadName && record.readName == traceReadName
            
            if(isTraceRead) {
                log.info "TRACE [$traceReadName]: Processing read at ${record.referenceName}:${record.alignmentStart}"
            }
            
            ReadWrapper wrapper = new ReadWrapper(
                record: record,
                sequenceNumber: seqNum++
            )
            
            // Fast filter - skip obvious non-candidates
            if (!hasSoftClips(record)) {
                if(isTraceRead) {
                    log.info "TRACE [$traceReadName]: REJECTED - No soft clips found"
                }
                wrapper.skipProcessing = true
                wrapper.processed = true
            }
            else if(!lengthsRoughlyMatch(record)) {
                if(isTraceRead) {
                    def info = detector.extractSoftClipInfo(record)
                    log.info "TRACE [$traceReadName]: REJECTED - Lengths don't match (soft clip: ${info.totalSoftClip}, aligned: ${info.alignedBases})"
                }
                wrapper.skipProcessing = true
                wrapper.processed = true
            }
            else {
                candidateReads.incrementAndGet()
                if(isTraceRead) {
                    log.info "TRACE [$traceReadName]: CANDIDATE - Passed fast filter, sending to worker for alignment check"
                }
            }
            
            // Add to both queues - workers read from input, writer reads from output
            inputQueue.put(wrapper)
            outputQueue.put(wrapper)
            
            long count = totalReads.incrementAndGet()
            
            // Print progress every PROGRESS_INTERVAL reads
            if(count % PROGRESS_INTERVAL == 0) {
                printProgress(count)
            }
        }
        
        if(region) {
            bam.withIterator(region, processRecord)
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
                
                if (wrapper.skipProcessing) {
                    // Already marked processed by fast filter - nothing to do
                    continue
                }
                
                boolean isTraceRead = traceReadName && wrapper.record.readName == traceReadName
                
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
                    wrapper.record = applyAction(wrapper.record)
                    
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
    SAMRecord applyAction(SAMRecord record) {
        String action = opts['action'] ?: 'reject'
        
        switch(action) {
            case 'reject':
                rejectedReads.incrementAndGet()
                return null  // Signal to skip this read
                
            case 'trim':
                trimmedReads.incrementAndGet()
                return trimSoftClips(record)
                
            case 'keep':
                return record  // Keep as-is
                
            default:
                log.warning "Unknown action: $action, defaulting to reject"
                rejectedReads.incrementAndGet()
                return null
        }
    }
    
    /**
     * Remove soft clips from a read, keeping only the aligned portion
     */
    @CompileStatic
    SAMRecord trimSoftClips(SAMRecord record) {
        def info = detector.extractSoftClipInfo(record)
        
        // Create a new CIGAR without soft clips
        def newCigar = record.cigar.cigarElements
            .findAll { it.operator != CigarOperator.SOFT_CLIP }
        
        if (newCigar.isEmpty()) {
            return null  // Nothing left after trimming
        }
        
        // Update the record
        record.cigar = new htsjdk.samtools.Cigar(newCigar)
        
        // Trim the bases and qualities
        int trimStart = info.leadingSoftClip
        int trimEnd = record.readLength - info.trailingSoftClip
        
        if (trimStart < trimEnd) {
            byte[] newBases = Arrays.copyOfRange(record.readBases, trimStart, trimEnd)
            byte[] newQuals = Arrays.copyOfRange(record.baseQualities, trimStart, trimEnd)
            
            record.readBases = newBases
            record.baseQualities = newQuals
        }
        
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
        
        log.info String.format(
            "Progress: %,d reads | %,d candidates (%.2f%%) | %,d duplex (%.2f%%) | %,d rejected | %,d trimmed",
            count, candidates, candidatePercent, duplex, duplexPercent, rejected, trimmed
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
}
