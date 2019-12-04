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

import cc.redberry.pipe.OutputPort;
import cc.redberry.pipe.OutputPortCloseable;
import cc.redberry.pipe.util.CountLimitingOutputPort;
import com.milaboratory.cli.PipelineConfiguration;
import com.milaboratory.mixcr.assembler.CloneAssemblerParameters;
import com.milaboratory.mixcr.vdjaligners.VDJCAlignerParameters;
import com.milaboratory.primitivio.PipeDataInputReader;
import com.milaboratory.primitivio.PrimitivI;
import com.milaboratory.primitivio.blocks.PrimitivIHybrid;
import com.milaboratory.util.CanReportProgress;
import io.repseq.core.GeneFeature;
import io.repseq.core.VDJCGene;
import io.repseq.core.VDJCLibraryRegistry;

import java.io.DataInput;
import java.io.DataInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Reader of CLNA file format.
 */
public final class ClnAReader extends PipelineConfigurationReaderMiXCR implements AutoCloseable {
    final PrimitivIHybrid input;

    // Index data

    final long firstClonePosition;
    // Index contain two additional records:
    //  - first = position of alignment block with cloneIndex == -1
    //  - last = position of the last alignments block end
    final long[] index;
    final long[] counts;
    final long totalAlignmentsCount;

    // From constructor

    final VDJCLibraryRegistry libraryRegistry;

    // Read form file header

    final PipelineConfiguration configuration;
    final VDJCAlignerParameters alignerParameters;
    final CloneAssemblerParameters assemblerParameters;

    final List<VDJCGene> genes;

    final int numberOfClones;

    // Meta data (also from header)

    final String versionInfo;

    public ClnAReader(Path path, VDJCLibraryRegistry libraryRegistry) throws IOException {
        this.input = new PrimitivIHybrid(ForkJoinPool.commonPool(), path);

        this.libraryRegistry = libraryRegistry;

        // File beginning
        String magicString;
        try (PrimitivI ii = this.input.beginPrimitivI()) {
            // Reading magic string
            byte[] magicBytes = new byte[ClnAWriter.MAGIC_LENGTH];
            ii.readFully(magicBytes);
            magicString = new String(magicBytes, StandardCharsets.US_ASCII);

            // Reading number of clones
            this.numberOfClones = ii.readInt();
        }

        // File ending
        long indexBegin;
        try (PrimitivI pi = this.input.beginRandomAccessPrimitivI(-IOUtil.END_MAGIC_LENGTH - 16)) {
            // Reading key file offsets from last 16 bytes of the file
            this.firstClonePosition = pi.readLong();
            indexBegin = pi.readLong();

            // Checking file consistency
            byte[] endMagic = new byte[IOUtil.END_MAGIC_LENGTH];
            pi.readFully(endMagic);
            if (!Arrays.equals(IOUtil.getEndMagicBytes(), endMagic))
                throw new RuntimeException("File is corrupted.");
        }

        // Step back
        try (PrimitivI pi = this.input.beginRandomAccessPrimitivI(indexBegin)) {
            // Reading index data
            this.index = new long[numberOfClones + 2];
            this.counts = new long[numberOfClones + 2];
            long previousValue = 0;
            long totalAlignmentsCount = 0L;
            for (int i = 0; i < numberOfClones + 2; i++) {
                previousValue = index[i] = previousValue + pi.readVarLong();
                totalAlignmentsCount += counts[i] = pi.readVarLong();
            }
            this.totalAlignmentsCount = totalAlignmentsCount;
        }

        // Returning to the file begin
        try (PrimitivI pi = this.input.beginPrimitivI(true)) {
            switch (magicString) {
                case ClnAWriter.MAGIC:
                    break;
                default:
                    throw new IllegalStateException();
            }

            this.versionInfo = pi.readUTF();
            this.configuration = pi.readObject(PipelineConfiguration.class);
            this.alignerParameters = pi.readObject(VDJCAlignerParameters.class);
            this.assemblerParameters = pi.readObject(CloneAssemblerParameters.class);
            this.genes = IOUtil.stdVDJCPrimitivIStateInit(pi, this.alignerParameters, libraryRegistry);
        }
    }

    public ClnAReader(String path, VDJCLibraryRegistry libraryRegistry) throws IOException {
        this(Paths.get(path), libraryRegistry);
    }

    @Override
    public PipelineConfiguration getPipelineConfiguration() {
        return configuration;
    }

    /**
     * Aligner parameters
     */
    public VDJCAlignerParameters getAlignerParameters() {
        return alignerParameters;
    }

    /**
     * Clone assembler parameters
     */
    public CloneAssemblerParameters getAssemblerParameters() {
        return assemblerParameters;
    }

    public GeneFeature[] getAssemblingFeatures() {
        return assemblerParameters.getAssemblingFeatures();
    }

    /**
     * Returns number of clones in the file
     */
    public int numberOfClones() {
        return numberOfClones;
    }

    public List<VDJCGene> getGenes() {
        return genes;
    }

    /**
     * Returns total number of alignments in the file, including unassembled.
     */
    public long numberOfAlignments() {
        return totalAlignmentsCount;
    }

    /**
     * Returns number of alignments contained in particular clone
     *
     * @param cloneIndex clone index
     * @return number of alignments
     */
    public long numberOfAlignmentsInClone(int cloneIndex) {
        return counts[cloneIndex + 1];
    }

    /**
     * MiXCR version this file was produced with.
     */
    public String getVersionInfo() {
        return versionInfo;
    }

    /**
     * Read clone set completely
     */
    public CloneSet readCloneSet() throws IOException {
        // TODO change to Blocks???

        // Reading clones
        int count = numberOfClones();
        List<Clone> clones = new ArrayList<>(count);

        try (PrimitivI pi = input.beginRandomAccessPrimitivI(firstClonePosition)) {
            for (int i = 0; i < count; i++)
                clones.add(pi.readObject(Clone.class));
        }

        return new CloneSet(clones, genes, alignerParameters, assemblerParameters);
    }

    /**
     * Constructs output port to read clones one by one as a stream
     */
    public OutputPortCloseable<Clone> readClones() throws IOException {
        PrimitivI input = inputState.createPrimitivI(new InputDataStream(firstClonePosition, index[0]));

        return new PipeDataInputReader<>(Clone.class, input, numberOfClones());
    }

    /**
     * Constructs output port to read alignments for a specific clone, or read unassembled alignments block
     *
     * @param cloneIndex index of clone; -1 to read unassembled alignments
     */
    public OutputPortCloseable<VDJCAlignments> readAlignmentsOfClone(int cloneIndex) {
        return new CountLimitingOutputPort<>(new BasicVDJCAlignmentReader(
                new AlignmentsIO.FileChannelBufferReader(channel, index[cloneIndex + 1], index[cloneIndex + 2]),
                inputState, false), counts[cloneIndex + 1]);
    }

    /**
     * Constructs output port to read all alignments form the file. Alignments are sorted by cloneIndex.
     */
    public OutputPortCloseable<VDJCAlignments> readAllAlignments() {
        return new CountLimitingOutputPort<>(new BasicVDJCAlignmentReader(
                new AlignmentsIO.FileChannelBufferReader(channel, index[0], index[index.length - 1]),
                inputState, false), totalAlignmentsCount);
    }

    /**
     * Constructs output port to read all alignments that are attached to a clone. Alignments are sorted by cloneIndex.
     */
    public OutputPortCloseable<VDJCAlignments> readAssembledAlignments() {
        return new CountLimitingOutputPort<>(new BasicVDJCAlignmentReader(
                new AlignmentsIO.FileChannelBufferReader(channel, index[1], index[index.length - 1]),
                inputState, false), totalAlignmentsCount - counts[0]);
    }

    /**
     * Constructs output port to read alignments that are not attached to any clone. Alignments are sorted by
     * cloneIndex.
     *
     * Returns: readAlignmentsOfClone(-1)
     */
    public OutputPortCloseable<VDJCAlignments> readNotAssembledAlignments() {
        return readAlignmentsOfClone(-1);
    }

    /**
     * Constructs output port of CloneAlignments objects, that allows to get synchronised view on clone and it's
     * corresponding alignments
     */
    public CloneAlignmentsPort clonesAndAlignments() throws IOException {
        return new CloneAlignmentsPort();
    }

    public final class CloneAlignmentsPort
            implements OutputPort<CloneAlignments>, CanReportProgress {
        private final AtomicLong processedAlignments = new AtomicLong();
        private final CloneSet fakeCloneSet;
        private final PipeDataInputReader<Clone> clones;
        volatile boolean isFinished = false;

        CloneAlignmentsPort() throws IOException {
            PrimitivI input = inputState.createPrimitivI(new InputDataStream(firstClonePosition, index[0]));
            this.clones = new PipeDataInputReader<>(Clone.class, input, numberOfClones);
            this.fakeCloneSet = new CloneSet(Collections.EMPTY_LIST, genes, alignerParameters, assemblerParameters);
        }

        @Override
        public CloneAlignments take() {
            Clone clone = clones.take();
            if (clone == null) {
                isFinished = true;
                return null;
            }
            clone.setParentCloneSet(fakeCloneSet);
            CloneAlignments result = new CloneAlignments(clone, clone.id);
            processedAlignments.addAndGet(result.alignmentsCount);
            return result;
        }

        @Override
        public double getProgress() {
            return 1.0 * processedAlignments.get() / totalAlignmentsCount;
        }

        @Override
        public boolean isFinished() {
            return isFinished;
        }
    }

    /**
     * Clone and alignments it was formed form
     */
    public final class CloneAlignments {
        /**
         * Clone
         */
        public final Clone clone;
        final int cloneId;
        final long alignmentsCount;

        CloneAlignments(Clone clone, int cloneId) {
            this.clone = clone;
            this.cloneId = cloneId;
            this.alignmentsCount = counts[cloneId + 1];
        }

        /**
         * Alignments
         */
        public OutputPort<VDJCAlignments> alignments() {
            return readAlignmentsOfClone(cloneId);
        }
    }

    @Override
    public void close() throws IOException {
        channel.close();
    }

    /**
     * FileChannel -> DataInput adapter that can be constructed for an arbitrary file position.
     *
     * Implemented using ByteBuffer.
     *
     * Thread-unsafe.
     */
    private class InputDataStream implements DataInput {
        private final long to;
        private final ByteBuffer buffer;
        private long lastPosition;

        InputDataStream(long from, long to) throws IOException {
            this.to = to;
            this.buffer = ByteBuffer.allocate(chunkSize);
            this.lastPosition = from;

            // Initially buffer is empty
            this.buffer.limit(0);

            // Filling first chunk of data
            if (from < to)
                fillBuffer();
        }

        void fillBuffer() throws IOException {
            // Number of bytes to read from file
            int size = (int) Math.min(chunkSize - buffer.remaining(), to - lastPosition);

            // Checking state
            if (size == 0)
                throw new IllegalArgumentException("No more bytes.");

            // Saving remaining bytes
            ByteBuffer remaining = buffer.slice();

            // Reset buffer state
            buffer.flip();

            // Setting new limit
            buffer.limit(size + remaining.limit());

            // Transferring remaining buffer bytes
            buffer.put(remaining);

            // Reading content form file
            int read = channel.read(buffer, lastPosition);

            // Flipping buffer
            buffer.flip();

            if (read != size)
                throw new IOException("Wrong block positions.");

            // Advancing last position
            this.lastPosition += read;
        }

        void ensureBuffer(int requiredSize) throws IOException {
            if (requiredSize > chunkSize)
                throw new IllegalArgumentException("Can't read this many bytes.");
            if (buffer.remaining() < requiredSize)
                fillBuffer();
        }

        @Override
        public void readFully(byte[] b) throws IOException {
            readFully(b, 0, b.length);
        }

        @Override
        public void readFully(byte[] b, int off, int len) throws IOException {
            do {
                int l = Math.min(chunkSize, len);
                ensureBuffer(l);
                buffer.get(b, off, l);
                off += l;
                len -= l;
            } while (len != 0);
        }

        @Override
        public int skipBytes(int n) throws IOException {
            ensureBuffer(n);
            buffer.position(buffer.position() + n);
            return n;
        }

        @Override
        public boolean readBoolean() throws IOException {
            byte b = buffer.get();
            if (b == 1)
                return true;
            else if (b == 0)
                return false;
            else
                throw new IOException("Illegal file format, can't deserialize boolean.");
        }

        @Override
        public byte readByte() throws IOException {
            ensureBuffer(1);
            return buffer.get();
        }

        @Override
        public int readUnsignedByte() throws IOException {
            ensureBuffer(1);
            return 0xFF & buffer.get();
        }

        @Override
        public short readShort() throws IOException {
            ensureBuffer(2);
            return (short) buffer.getChar();
        }

        @Override
        public int readUnsignedShort() throws IOException {
            ensureBuffer(2);
            return 0xFFFF & buffer.getChar();
        }

        @Override
        public char readChar() throws IOException {
            ensureBuffer(2);
            return buffer.getChar();
        }

        @Override
        public int readInt() throws IOException {
            ensureBuffer(4);
            return buffer.getInt();
        }

        @Override
        public long readLong() throws IOException {
            ensureBuffer(8);
            return buffer.getLong();
        }

        @Override
        public float readFloat() throws IOException {
            ensureBuffer(4);
            return buffer.getFloat();
        }

        @Override
        public double readDouble() throws IOException {
            ensureBuffer(8);
            return buffer.getDouble();
        }

        @Override
        public String readLine() throws IOException {
            throw new UnsupportedOperationException();
        }

        @Override
        public String readUTF() throws IOException {
            return DataInputStream.readUTF(this);
        }
    }
}
