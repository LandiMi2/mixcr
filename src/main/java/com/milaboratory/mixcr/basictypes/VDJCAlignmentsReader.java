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

import cc.redberry.pipe.OutputPortCloseable;
import com.milaboratory.cli.PipelineConfiguration;
import com.milaboratory.mixcr.vdjaligners.VDJCAlignerParameters;
import com.milaboratory.primitivio.PrimitivI;
import com.milaboratory.primitivio.blocks.PrimitivIBlocks;
import com.milaboratory.primitivio.blocks.PrimitivIBlocksStats;
import com.milaboratory.primitivio.blocks.PrimitivIHybrid;
import com.milaboratory.util.CanReportProgress;
import io.repseq.core.GeneFeature;
import io.repseq.core.GeneType;
import io.repseq.core.VDJCGene;
import io.repseq.core.VDJCLibraryRegistry;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ForkJoinPool;

import static com.milaboratory.mixcr.basictypes.VDJCAlignmentsWriter.*;

public final class VDJCAlignmentsReader extends PipelineConfigurationReaderMiXCR implements
        OutputPortCloseable<VDJCAlignments>,
        CanReportProgress {
    public static final int DEFAULT_CONCURRENCY = 4;
    public static final int DEFAULT_READ_AHEAD_BLOCKS = 5;

    final PrimitivIHybrid input;
    final int readAheadBlocks, concurrency;
    final VDJCLibraryRegistry vdjcRegistry;

    PrimitivIBlocks<VDJCAlignments>.Reader reader;

    VDJCAlignerParameters parameters;
    PipelineConfiguration pipelineConfiguration;
    List<VDJCGene> usedGenes;

    String versionInfo;
    String magic;
    long counter = 0;
    long numberOfReads = -1;
    boolean closed = false;
    final long size;

    public VDJCAlignmentsReader(String fileName) throws IOException {
        this(fileName, VDJCLibraryRegistry.getDefault());
    }

    public VDJCAlignmentsReader(String fileName, VDJCLibraryRegistry vdjcRegistry) throws IOException {
        this(fileName, vdjcRegistry, DEFAULT_CONCURRENCY);
    }

    public VDJCAlignmentsReader(String fileName, VDJCLibraryRegistry vdjcRegistry, int concurrency) throws IOException {
        this(Paths.get(fileName), vdjcRegistry, concurrency);
    }

    public VDJCAlignmentsReader(File file) throws IOException {
        this(file, VDJCLibraryRegistry.getDefault());
    }

    public VDJCAlignmentsReader(File file, VDJCLibraryRegistry vdjcRegistry) throws IOException {
        this(file, vdjcRegistry, DEFAULT_CONCURRENCY);
    }

    public VDJCAlignmentsReader(File file, VDJCLibraryRegistry vdjcRegistry, int concurrency) throws IOException {
        this(file.toPath(), vdjcRegistry, concurrency);
    }

    public VDJCAlignmentsReader(Path path) throws IOException {
        this(path, VDJCLibraryRegistry.getDefault());
    }

    public VDJCAlignmentsReader(Path path, VDJCLibraryRegistry vdjcRegistry) throws IOException {
        this(path, vdjcRegistry, DEFAULT_CONCURRENCY);
    }

    public VDJCAlignmentsReader(Path path, VDJCLibraryRegistry vdjcRegistry, int concurrency) throws IOException {
        this(path, vdjcRegistry, concurrency, ForkJoinPool.commonPool());
    }

    public VDJCAlignmentsReader(Path path, VDJCLibraryRegistry vdjcRegistry, int concurrency, ExecutorService executor) throws IOException {
        this(path, vdjcRegistry, concurrency, executor, DEFAULT_READ_AHEAD_BLOCKS);
    }

    public VDJCAlignmentsReader(Path path, VDJCLibraryRegistry vdjcRegistry, int concurrency, ExecutorService executor, int readAheadBlocks) throws IOException {
        this.input = new PrimitivIHybrid(executor, path);
        this.readAheadBlocks = readAheadBlocks;
        this.vdjcRegistry = vdjcRegistry;
        this.concurrency = concurrency;
        this.size = Files.size(path);
    }

    // public void init() {
    //     init(null);
    // }

    public void init() {
        if (reader != null)
            return;

        try (final PrimitivI i = input.beginPrimitivI(true)) {
            assert MAGIC_BYTES.length == MAGIC_LENGTH;
            byte[] magic = new byte[MAGIC_LENGTH];
            i.readFully(magic);
            String magicString = new String(magic);
            this.magic = magicString;

            // SerializersManager serializersManager = input.getSerializersManager();
            switch (magicString) {
                case MAGIC:
                    break;
                default:
                    throw new RuntimeException("Unsupported file format; .vdjca file of version " + new String(magic)
                            + " while you are running MiXCR " + MAGIC);
            }

            versionInfo = i.readUTF();

            parameters = i.readObject(VDJCAlignerParameters.class);
            pipelineConfiguration = i.readObject(PipelineConfiguration.class);

            this.usedGenes = IOUtil.stdVDJCPrimitivIStateInit(i, parameters, vdjcRegistry);
        }

        this.reader = input.beginPrimitivIBlocks(VDJCAlignments.class, concurrency, readAheadBlocks);
    }

    public PrimitivIBlocksStats getStats() {
        if (reader == null)
            return null;
        return reader.getParent().getStats();
    }

    public synchronized VDJCAlignerParameters getParameters() {
        init();
        return parameters;
    }

    public synchronized List<VDJCGene> getUsedGenes() {
        init();
        return usedGenes;
    }

    @Override
    public synchronized PipelineConfiguration getPipelineConfiguration() {
        init();
        return pipelineConfiguration;
    }

    /**
     * Returns information about version of MiXCR which produced this file.
     *
     * @return information about version of MiXCR which produced this file
     */
    public String getVersionInfo() {
        init();
        return versionInfo;
    }

    /**
     * Returns magic bytes of this file.
     *
     * @return magic bytes of this file
     */
    public String getMagic() {
        init();
        return magic;
    }

    public long getNumberOfReads() {
        return numberOfReads;
    }

    @Override
    public double getProgress() {
        if (size == 0)
            return Double.NaN;
        return (1.0 * input.getPosition()) / size;
    }

    @Override
    public boolean isFinished() {
        return closed || input.getPosition() == size;
    }

    @Override
    public synchronized void close() {
        close(false);
    }

    private void close(boolean onEnd) {
        if (closed)
            return;

        try {
            // Closing blocked reader
            reader.close();

            // If all alignments are read
            // footer with number of reads processed to produce this
            // file can be read form the stream.
            if (onEnd)
                try (final PrimitivI i = input.beginPrimitivI()) {
                    numberOfReads = i.readLong();
                }

            input.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            closed = true;
        }
    }

    @Override
    public synchronized VDJCAlignments take() {
        if (closed)
            return null;

        init();

        VDJCAlignments al = reader.take();

        if (al == null) {
            close(true);
            return null;
        }

        return al.setAlignmentsIndex(counter++);
    }

    // /**
    //  * Produce reader that uses the same reference for geneFeatures.
    //  *
    //  * @param reader     target reader
    //  * @param parameters parameters to take reference from
    //  */
    // public static void initGeneFeatureReferencesFrom(VDJCAlignmentsReader reader, VDJCAlignerParameters parameters) {
    //     Map<GeneFeature, GeneFeature> featureRefs = new HashMap<>();
    //     for (GeneType gt : GeneType.VDJC_REFERENCE) {
    //         GeneFeature f = parameters.getFeatureToAlign(gt);
    //         featureRefs.put(f, f);
    //     }
    //     reader.init(featureRefs);
    // }
    //
    // /**
    //  * Produce reader that uses the same reference for geneFeatures.
    //  *
    //  * @param reader       target reader
    //  * @param sourceReader reader to take reference from
    //  */
    // public static void initGeneFeatureReferencesFrom(VDJCAlignmentsReader reader, VDJCAlignmentsReader sourceReader) {
    //     initGeneFeatureReferencesFrom(reader, sourceReader.getParameters());
    // }
}
