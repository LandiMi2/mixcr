/*
 * Copyright (c) 2014-2019, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact MiLaboratory LLC, which owns exclusive
 * rights for distribution of this program for commercial purposes, using the
 * following email address: licensing@milaboratory.com.
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */
package com.milaboratory.mixcr.basictypes;

import com.milaboratory.cli.AppVersionInfo;
import com.milaboratory.cli.PipelineConfiguration;
import com.milaboratory.mixcr.util.MiXCRVersionInfo;
import com.milaboratory.mixcr.vdjaligners.VDJCAligner;
import com.milaboratory.mixcr.vdjaligners.VDJCAlignerParameters;
import com.milaboratory.primitivio.PrimitivO;
import com.milaboratory.primitivio.blocks.PrimitivOBlocks;
import com.milaboratory.primitivio.blocks.PrimitivOBlocksStats;
import com.milaboratory.primitivio.blocks.PrimitivOHybrid;
import com.milaboratory.util.io.HasPosition;
import io.repseq.core.VDJCGene;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ForkJoinPool;

import static com.milaboratory.mixcr.basictypes.AlignmentsIO.DEFAULT_ALIGNMENTS_IN_BLOCK;

public final class VDJCAlignmentsWriter implements VDJCAlignmentsWriterI, HasPosition {
    public static final int DEFAULT_ENCODER_THREADS = 3;
    static final String MAGIC_V15 = "MiXCR.VDJC.V15";
    static final String MAGIC = MAGIC_V15;
    static final int MAGIC_LENGTH = 14;
    static final byte[] MAGIC_BYTES = MAGIC.getBytes(StandardCharsets.US_ASCII);

    /**
     * This number will be added to the end of the file to report number of processed read to the following processing
     * steps. Mainly for informative statistics reporting.
     */
    long numberOfProcessedReads = -1;

    /** Writer settings */
    final int encoderThreads, alignmentsInBlock;

    /** Main block output class */
    final PrimitivOHybrid output;

    /** Initialized after headers */
    volatile PrimitivOBlocks<VDJCAlignments>.Writer writer;

    boolean closed = false;

    public VDJCAlignmentsWriter(String fileName) throws IOException {
        this(fileName, DEFAULT_ENCODER_THREADS, DEFAULT_ALIGNMENTS_IN_BLOCK);
    }

    public VDJCAlignmentsWriter(String fileName, int encoderThreads, int alignmentsInBlock) throws IOException {
        this(fileName, encoderThreads, alignmentsInBlock, ForkJoinPool.commonPool());
    }

    public VDJCAlignmentsWriter(String fileName, int encoderThreads, int alignmentsInBlock, ExecutorService executor) throws IOException {
        this(new PrimitivOHybrid(executor, Paths.get(fileName)), encoderThreads, alignmentsInBlock);
    }

    public VDJCAlignmentsWriter(File file) throws IOException {
        this(file, DEFAULT_ENCODER_THREADS, DEFAULT_ALIGNMENTS_IN_BLOCK);
    }

    public VDJCAlignmentsWriter(File file, int encoderThreads, int alignmentsInBlock) throws IOException {
        this(file, encoderThreads, alignmentsInBlock, ForkJoinPool.commonPool());
    }

    public VDJCAlignmentsWriter(File file, int encoderThreads, int alignmentsInBlock, ExecutorService executor) throws IOException {
        this(new PrimitivOHybrid(executor, file.toPath()), encoderThreads, alignmentsInBlock);
    }

    public VDJCAlignmentsWriter(PrimitivOHybrid output, int encoderThreads, int alignmentsInBlock) {
        this.output = output;
        this.encoderThreads = encoderThreads;
        this.alignmentsInBlock = alignmentsInBlock;
    }

    @Override
    public void setNumberOfProcessedReads(long numberOfProcessedReads) {
        this.numberOfProcessedReads = numberOfProcessedReads;
    }

    public void header(VDJCAlignmentsReader reader, PipelineConfiguration pipelineConfiguration) {
        header(reader.getParameters(), reader.getUsedGenes(), pipelineConfiguration);
    }

    public void header(VDJCAligner aligner, PipelineConfiguration pipelineConfiguration) {
        header(aligner.getParameters(), aligner.getUsedGenes(), pipelineConfiguration);
    }

    /** History to write in the header */
    private PipelineConfiguration pipelineConfiguration = null;

    public synchronized void setPipelineConfiguration(PipelineConfiguration configuration) {
        if (pipelineConfiguration == null)
            pipelineConfiguration = configuration;
        else if (!configuration.equals(this.pipelineConfiguration))
            throw new IllegalStateException();
    }

    @Override
    public void header(VDJCAlignerParameters parameters, List<VDJCGene> genes,
                       PipelineConfiguration ppConfiguration) {
        if (parameters == null || genes == null)
            throw new IllegalArgumentException();

        if (writer != null)
            throw new IllegalStateException("Header already written.");

        // Writing meta data using raw stream for easy reconstruction with simple tools like hex viewers
        try (PrimitivO o = output.beginPrimitivO(true)) {
            // Writing magic bytes
            assert MAGIC_BYTES.length == MAGIC_LENGTH;
            o.write(MAGIC_BYTES);

            // Writing version information
            o.writeUTF(
                    MiXCRVersionInfo.get().getVersionString(
                            AppVersionInfo.OutputType.ToFile));

            // Writing parameters
            o.writeObject(parameters);

            // Writing history
            if (ppConfiguration != null)
                this.pipelineConfiguration = ppConfiguration;
            o.writeObject(pipelineConfiguration);

            IOUtil.stdVDJCPrimitivOStateInit(o, genes, parameters);
        }

        writer = output.beginPrimitivOBlocks(encoderThreads, alignmentsInBlock);
    }

    @Override
    public long getPosition() {
        return output.getPosition();
    }

    public int getEncodersCount() {
        if (writer == null)
            return 0;
        return writer.getParent().getStats().getConcurrency();
    }

    public int getBusyEncoders() {
        if (writer == null)
            return 0;
        PrimitivOBlocksStats stats = writer.getParent().getStats();
        return stats.getOngoingSerdes() + stats.getOngoingIOOps();
    }

    @Override
    public synchronized void write(VDJCAlignments alignment) {
        if (writer == null)
            throw new IllegalStateException("Header not initialized.");

        if (alignment == null)
            throw new NullPointerException();

        writer.write(alignment);
    }

    @Override
    public synchronized void close() {
        try {
            if (!closed) {
                writer.close(); // This will also write stream termination symbol/block to the stream

                try (PrimitivO o = output.beginPrimitivO()) {
                    // [ numberOfProcessedReads : long ]
                    byte[] footer = new byte[8];

                    // Number of processed reads is known only in the end of analysis
                    // Writing it as last piece of information in the stream
                    AlignmentsIO.writeLongBE(numberOfProcessedReads, footer, 0);
                    o.write(footer);

                    // Writing end-magic as a file integrity sign
                    o.write(IOUtil.getEndMagicBytes());
                } finally {
                    closed = true;
                    output.close();
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public boolean isClosed() {
        return closed;
    }
}
