package gngs.tools

import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.atomic.AtomicLong
import java.util.concurrent.ConcurrentHashMap

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.samtools.CigarElement
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
    AtomicLong removedSecondaries = new AtomicLong(0)
    
    // Trace logging
    String traceReadName = null
    
    // Secondary alignment tracking
    static final int MAX_SECONDARY_BUFFER = 1000
    
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
                sequenceNumber: seqNum++,
                originalReadName: record.readName
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
        
        return detector.lengthsMatch(info.leadingSoftClip, info.alignedBases) ||
               detector.lengthsMatch(info.trailingSoftClip, info.alignedBases)
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
                    boolean leadingMatches = threadDetector.lengthsMatch(duplexInfo.leadingSoftClip, duplexInfo.alignedBases)
                    boolean trailingMatches = threadDetector.lengthsMatch(duplexInfo.trailingSoftClip, duplexInfo.alignedBases)
                    if (trailingMatches && leadingMatches) {
                        wrapper.duplexOnTrailingEnd = (duplexInfo.trailingSoftClip >= duplexInfo.leadingSoftClip)
                    } else {
                        wrapper.duplexOnTrailingEnd = trailingMatches
                    }
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
     * Start the output thread that writes processed reads in order from the output queue.
     * 
     * When removeSecondaries is enabled, this thread tracks which read names have been
     * identified as duplex. If a secondary/supplementary alignment is encountered whose
     * primary hasn't been resolved yet, the output stream is buffered until the primary
     * is found (which should be very nearby in a coordinate-sorted BAM since the secondary
     * aligns to approximately the same position as the primary).
     */
    Thread startOutputThread(LinkedBlockingQueue<ReadWrapper> outputQueue, 
                            ReadWrapper POISON, String outputFile, SAM bam) {
        
        boolean removeSecondaries = opts['removeSecondaries'] ? true : false
        
        Thread.start {
            log.info "Writing output to: $outputFile"
            if (removeSecondaries) {
                log.info "Secondary alignment removal enabled"
            }
            
            bam.withWriter(outputFile, true) { SAMFileWriter writer ->
                
                // Track resolved read names
                Set<String> duplexReadNames = new HashSet<>()
                Set<String> resolvedNonDuplex = new HashSet<>()
                
                // Eviction tracking - remove old entries to bound memory
                long evictionCounter = 0
                final long EVICTION_INTERVAL = 50000
                
                while (true) {
                    ReadWrapper wrapper = outputQueue.take()
                    
                    if (wrapper.sequenceNumber == -1) {
                        break
                    }
                    
                    // Wait until this specific record is processed
                    while (!wrapper.processed) {
                        Thread.sleep(1)
                    }
                    
                    // null record means it was rejected (but we still need to track the name)
                    if (wrapper.record == null) {
                        if (removeSecondaries && wrapper.isDuplex) {
                            duplexReadNames.add(wrapper.originalReadName)
                        }
                        continue
                    }
                    
                    SAMRecord record = wrapper.record
                    
                    if (removeSecondaries && record.isSecondaryOrSupplementary()) {
                        String readName = record.readName
                        if (duplexReadNames.contains(readName)) {
                            // Primary was already resolved as duplex - discard
                            removedSecondaries.incrementAndGet()
                            continue
                        } else if (resolvedNonDuplex.contains(readName)) {
                            // Primary was already resolved as non-duplex - write
                            writer.addAlignment(record)
                        } else {
                            // Primary not yet seen - buffer until resolved
                            List<ReadWrapper> localBuffer = new ArrayList<>()
                            localBuffer.add(wrapper)
                            boolean resolved = false
                            
                            while (!resolved) {
                                ReadWrapper next = outputQueue.take()
                                
                                if (next.sequenceNumber == -1) {
                                    // End of stream - flush buffer and exit
                                    for (ReadWrapper buffered : localBuffer) {
                                        if (buffered.record != null) {
                                            writer.addAlignment(buffered.record)
                                        }
                                    }
                                    resolved = true
                                    // Re-signal poison for the outer loop
                                    outputQueue.put(next)
                                    break
                                }
                                
                                while (!next.processed) {
                                    Thread.sleep(1)
                                }
                                
                                localBuffer.add(next)
                                
                                // Check if this is the primary that resolves our secondary
                                if (next.record != null && !next.record.isSecondaryOrSupplementary()
                                    && next.record.readName == readName) {
                                    if (next.isDuplex) {
                                        duplexReadNames.add(readName)
                                    } else {
                                        resolvedNonDuplex.add(readName)
                                    }
                                    resolved = true
                                } else if (next.record == null && next.isDuplex 
                                           && next.originalReadName == readName) {
                                    // Primary was rejected (record set to null)
                                    duplexReadNames.add(readName)
                                    resolved = true
                                }
                                
                                // Safety valve: if buffer grows too large, assume non-duplex
                                if (!resolved && localBuffer.size() > MAX_SECONDARY_BUFFER) {
                                    log.warning "Secondary buffer exceeded ${MAX_SECONDARY_BUFFER} for read ${readName}, assuming non-duplex"
                                    resolvedNonDuplex.add(readName)
                                    resolved = true
                                }
                            }
                            
                            // Flush local buffer in order
                            for (ReadWrapper buffered : localBuffer) {
                                if (buffered.record == null) {
                                    continue
                                }
                                if (buffered.record.isSecondaryOrSupplementary()
                                    && duplexReadNames.contains(buffered.record.readName)) {
                                    removedSecondaries.incrementAndGet()
                                    continue
                                }
                                writer.addAlignment(buffered.record)
                            }
                        }
                    } else {
                        // Primary read or removeSecondaries disabled
                        if (removeSecondaries) {
                            if (wrapper.isDuplex) {
                                duplexReadNames.add(record.readName)
                            } else {
                                resolvedNonDuplex.add(record.readName)
                            }
                        }
                        writer.addAlignment(record)
                    }
                    
                    // Periodic eviction of old resolved names to bound memory
                    if (removeSecondaries && (++evictionCounter % EVICTION_INTERVAL == 0)) {
                        // Keep duplexReadNames (small set, needed for distant secondaries)
                        // Clear resolvedNonDuplex (large set, distant secondaries are safe to write)
                        resolvedNonDuplex.clear()
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
        
        List<CigarElement> cigarElements = record.cigar.cigarElements
        List<CigarElement> newCigar = new ArrayList<CigarElement>()
        int trimStart = 0
        int trimEnd = record.readLength
        
        if (duplexOnTrailingEnd) {
            // Remove trailing soft clip (and any trailing hard clip) only
            trimEnd = record.readLength - info.trailingSoftClip
            
            // Build new CIGAR: keep everything except trailing S and H
            int lastNonTrailingIdx = cigarElements.size() - 1
            // Walk backwards past trailing H and S elements
            while (lastNonTrailingIdx >= 0) {
                CigarOperator op = cigarElements[lastNonTrailingIdx].operator
                if (op == CigarOperator.SOFT_CLIP || op == CigarOperator.HARD_CLIP) {
                    lastNonTrailingIdx--
                } else {
                    break
                }
            }
            // Copy elements up to and including the last non-trailing element
            for (int i = 0; i <= lastNonTrailingIdx; i++) {
                newCigar.add(cigarElements[i])
            }
        } else {
            // Remove leading soft clip (and any leading hard clip) only
            trimStart = info.leadingSoftClip
            
            // Build new CIGAR: keep everything except leading S and H
            int firstNonLeadingIdx = 0
            // Walk forward past leading H and S elements
            while (firstNonLeadingIdx < cigarElements.size()) {
                CigarOperator op = cigarElements[firstNonLeadingIdx].operator
                if (op == CigarOperator.SOFT_CLIP || op == CigarOperator.HARD_CLIP) {
                    firstNonLeadingIdx++
                } else {
                    break
                }
            }
            // Copy elements from the first non-leading element onwards
            for (int i = firstNonLeadingIdx; i < cigarElements.size(); i++) {
                newCigar.add(cigarElements[i])
            }
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
        
        long secondaries = removedSecondaries.get()
        
        log.info String.format(
            "Progress: %,d reads at %s | Queue: in=%d/%d (%.1f%%) out=%d/%d (%.1f%%) | %,d candidates (%.2f%%) | %,d duplex (%.2f%%) | %,d rejected | %,d trimmed | %,d sec. removed",
            count, position, inputQueueSize, totalCapacity, inputUtilization, outputQueueSize, totalCapacity, outputUtilization,
            candidates, candidatePercent, duplex, duplexPercent, rejected, trimmed, secondaries
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
        System.err.println "Secondaries removed:      ${removedSecondaries.get()}"
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
            removeSecondaries 'Remove secondary/supplementary alignments for duplex reads', required: false
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
    String originalReadName
    long sequenceNumber
    volatile boolean processed = false
    volatile boolean isDuplex = false
    volatile boolean skipProcessing = false
    volatile boolean duplexOnTrailingEnd = false
}
