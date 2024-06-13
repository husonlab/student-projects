package org.husonlab.fmhdist.util.experimental;

import jloda.util.FileUtils;
import org.husonlab.fmhdist.util.KMerCoordinates;
import org.husonlab.fmhdist.util.KMerIterator;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Iterator to extract all valid k-mers from a fasta file (raw, zipped or
 * gzipped), local file system or ftp url.
 * <p>
 * The iterator can skip ambiguous bases (N/n) and calculate the reverse
 * complement for the given k-mer.
 * <p>
 * With each new record in the file, a new k-mer is created, i.e. this iterator
 * won't return k-mers spanning multiple records.
 */
public class FastKMerIterator implements KMerIterator {
	private static boolean[] isLineContainingSkippableChar = new boolean[128];

	static {
		isLineContainingSkippableChar['\t'] = true;
		isLineContainingSkippableChar['>'] = true;
		isLineContainingSkippableChar['\n'] = true;
		isLineContainingSkippableChar['\r'] = true;
		isLineContainingSkippableChar[' '] = true;
	}

	private static final byte[] complementTable = new byte[128];

	static {
		complementTable['A'] = 'T';
		complementTable['T'] = 'A';
		complementTable['G'] = 'C';
		complementTable['C'] = 'G';
		complementTable['a'] = 't';
		complementTable['t'] = 'a';
		complementTable['g'] = 'c';
		complementTable['c'] = 'g';
	}

	private static final byte[] toUpperTable = new byte[128];

	static {
		toUpperTable['A'] = 'A';
		toUpperTable['T'] = 'T';
		toUpperTable['G'] = 'G';
		toUpperTable['C'] = 'C';
		toUpperTable['a'] = 'A';
		toUpperTable['t'] = 'T';
		toUpperTable['g'] = 'G';
		toUpperTable['c'] = 'C';
		toUpperTable['n'] = 'N';
		toUpperTable['N'] = 'N';
	}

	private final boolean[] isAmbiguousChar = new boolean[128];
	private final byte[] kmer;
	private final byte[] complement;
	private final byte[] preloaded_kmer;
	private final byte[] preloaded_complement;
	private final InputStream reader;
	private final int k;

	private byte nextByte;
	private boolean hasPreloaded;

	// Indices to keep track of the origin of the _current_ k-mer. The values
	// will always be copied from the preloaded variants (see below). The values
	// will be fed into the KMerCoordinates if the corresponding method
	// getCoordinates() is called. The separation enables to prepare the next
	// k-mer during a call to "next()" while still be able to load details of
	// the current kmer that is returned by the same call to next().
	private int recordIndexInFile = 0;
	private int skippedKmersInFile = 0;
	private int skippedKmersInRecord = 0;
	private int sequenceIndexInRecord = 0;
	private int sequenceIndexInFile = 0;

	// Indices to keep track of the origin of the _next_ k-mer. Those values
	// will actually be incremented as the stream is processed.Those will be
	// written to the current indices during the next call to next().
	private int preloadedRecordIndexInFile = 0;
	private int preloadedSkippedKmersInFile = 0;
	private int preloadedSkippedKmersInRecord = 0;
	private int preloadedSequenceIndexInRecord = 0;
	private int preloadedSequenceIndexInFile = 0;

	/**
	 * Create a new Iterator to extract the kmers from the given file.
	 *
	 * @param k        the size of the k-mer
	 * @param fileName the path to the file to fetch
	 * @param skipN    indicate if k-mers containing the letter "N" or "n" should
	 *                 be skipped. This is typically the case if the base is ambiguous for
	 *                 genomic sequences.
	 * @throws IOException
	 */
	public FastKMerIterator(int k, String fileName, boolean skipN) throws IOException {
		this(k, new BufferedInputStream(FileUtils.getInputStreamPossiblyZIPorGZIP(fileName)), skipN);
	}

	/**
	 * Create a new Iterator to extract the kmers from the given file.
	 *
	 * @param k      the size of the k-mer input strea
	 * @param reader the InputStreamReader to read the sequence from.
	 * @param skipN  indicate if k-mers containing the letter "N" or "n" should
	 *               be skipped. This is typically the case if the base is ambiguous for
	 *               genomic sequences.
	 * @throws IOException
	 */
	public FastKMerIterator(int k, InputStream reader, boolean skipN) throws IOException {
		this.k = k;
		this.kmer = new byte[k];
		this.preloaded_kmer = new byte[k];
		this.complement = new byte[k];
		this.preloaded_complement = new byte[k];

		this.isAmbiguousChar['N'] = skipN;
		this.isAmbiguousChar['n'] = skipN;

		this.reader = reader;
		this.nextByte = (byte) reader.read();
		this.preloadedRecordIndexInFile = -1;
		this.preload();
	}

	/**
	 * Returns the k Parameter of the Iterator.
	 *
	 * @return
	 */
	public int getK() {
		return this.k;
	}

	@Override
	public void close() throws IOException {
		this.reader.close();
	}

	@Override
	public boolean hasNext() {
		return this.nextByte != -1;
	}

	/**
	 * Returns true if the byte at this.nextByte is a character that is not
	 * known to be a Non-Sequence character.
	 *
	 * @return
	 */
	private boolean isSequenceChar() {
		return !isLineContainingSkippableChar[this.nextByte];
	}

	/**
	 * Returns true if the byte at this.nextByte is ">" which indicates the
	 * start of a new header (but only if this is also preceeded by "\n")
	 *
	 * @return
	 */
	private boolean isHeaderStart() {
		return this.nextByte == '>';
	}

	/**
	 * Proceeds to read bytes from the input stream until the start of the new
	 * line is reached. Only the preloadedByteCounter is incremented here.
	 *
	 * @throws IOException
	 */
	private void skipToNextLine() throws IOException {
		while (this.hasNext() && this.nextByte != '\n') {
			this.nextByte = (byte) this.reader.read();
		}
		if (this.hasNext()) {
			this.nextByte = (byte) this.reader.read();
		}
	}

	/**
	 * Whenever a new seqeunce header is reached in a FASTA file, we perform
	 * three steps:
	 * 1. Incremend the record counter (i.e. we are in the next sequence in the
	 * file)
	 * 2. Reset all Record-related counter to 0 (i.e. the next k-mer is the
	 * first in the record)
	 * 3. Skip to the next line.
	 *
	 * @throws IOException
	 */
	private void handleSequenceStart() throws IOException {
		this.preloadedRecordIndexInFile++;
		this.preloadedSequenceIndexInRecord = 0;
		this.preloadedSkippedKmersInRecord = 0;
		skipToNextLine();
	}

	/**
	 * Assuming that "this.nextByte" has NOT been handled yet, this function
	 * reads as many bytes as needed to fill in the first k-1 bytes of a NEW
	 * k-mer. This is needed if
	 * 1. A new sequence starts
	 * 2. An existing k-mer cannot be extended because the next character is
	 * ambiguous. In this case, the caller MUST already have loaded the byte
	 * AFTER the ambiguous nucleotide.
	 * <p>
	 * If this function needs to skip k-mers because it encounters ambiguous
	 * nucleotides, it updates the according counters.
	 *
	 * @throws IOException
	 */
	private void preload() throws IOException {
		this.hasPreloaded = true;
		int i = 0;
		while (this.hasNext()) {
			while (i < this.k - 1 && this.hasNext()) {
				if (isSequenceChar()) {
					if (isAmbiguousChar[this.nextByte]) {
						// we need to discard the previous i k-mers (0-based, thus
						// +1)
						this.preloadedSkippedKmersInFile += i + 1;
						this.preloadedSkippedKmersInRecord += i + 1;
						i = 0;
						this.nextByte = (byte) this.reader.read();
						break;
					}
					this.preloaded_kmer[i] = toUpperTable[this.nextByte];
					this.preloaded_complement[this.k - i - 1] = toUpperTable[complementTable[this.nextByte]];
					this.nextByte = (byte) this.reader.read();
					i++;
				} else {
					if (isHeaderStart()) {
						handleSequenceStart();
						i = 0;
					} else {
						skipToNextLine();
					}
				}
			}

			// Now, we need only to ensure that the last character is valid
			// If it is, return.
			// If it is a line break, read new character.
			if (i >= this.k - 1 && this.hasNext()) {
				if (this.isSequenceChar()) {
					if (!isAmbiguousChar[this.nextByte]) {
						return;
					} else {
						// Start all over again :(
						this.preloadedSkippedKmersInFile += i + 1;
						this.preloadedSkippedKmersInRecord += i + 1;
						i = 0;
						this.nextByte = (byte) this.reader.read();
						continue;
					}
				} else {
					if (this.isHeaderStart()) {
						this.handleSequenceStart();
						i = 0;
						continue;
					} else {
						this.skipToNextLine();
					}
				}
			}
		}
	}

	/**
	 * This updates all indices such that they equal the preloaded ones.
	 */
	private void copyIndices() {
		this.recordIndexInFile = this.preloadedRecordIndexInFile;
		this.sequenceIndexInFile = this.preloadedSequenceIndexInFile;
		this.sequenceIndexInRecord = this.preloadedSequenceIndexInRecord;
		this.skippedKmersInFile = this.preloadedSkippedKmersInFile;
		this.skippedKmersInRecord = this.preloadedSkippedKmersInRecord;
	}

	/**
	 * Returns the next k-mer and prepares a new one, if applicable. Only call
	 * this if "hasNext()" returns true.
	 * <p>
	 * After next() is called, calls to getReverseComplement() and
	 * getCoordinates() will return objects that belong to the k-mer that next()
	 * returned.
	 */
	@Override
	public byte[] next() {
		// First, finalize current k-mer
		if (this.hasPreloaded) {
			System.arraycopy(this.preloaded_kmer, 0, this.kmer, 0, this.k - 1);
			System.arraycopy(this.preloaded_complement, 1, this.complement, 1, this.k - 1);
			this.hasPreloaded = false;
		} else {
			System.arraycopy(this.kmer, 1, this.kmer, 0, this.k - 1);
			System.arraycopy(this.complement, 0, this.complement, 1, this.k - 1);
		}
		this.copyIndices();
		this.kmer[this.k - 1] = toUpperTable[this.nextByte];
		this.complement[0] = toUpperTable[complementTable[this.nextByte]];

		// Now, k-mer is finished. Time to prepare the next one!
		try {
			this.nextByte = (byte) reader.read();
			this.preloadedSequenceIndexInFile++;
			this.preloadedSequenceIndexInRecord++;

			boolean forceNextIteration = true;
			boolean needsPreload = false;

			// There are four possible cases to consider:
			// 1. The next character is a valid sequence character
			// 2. The next character is a ambiguous sequence character
			// 3. The next character is the start of a new header
			// 4. The next character is a \n
			while (hasNext() && forceNextIteration) {
				forceNextIteration = false;
				if (isSequenceChar()) {
					if (isAmbiguousChar[this.nextByte]) {
						// We need to skip the next k k-mers
						this.preloadedSkippedKmersInFile += this.k;
						this.preloadedSkippedKmersInRecord += this.k;
						this.nextByte = (byte) reader.read();
						needsPreload = true;
					}
				} else {
					if (isHeaderStart()) {
						handleSequenceStart();
						needsPreload = true;
					} else {
						skipToNextLine();
						forceNextIteration = true;
					}
				}

				if (needsPreload) {
					this.preload();
					// We do not need to break because in above, it is
					// impossible to needPreload AND forceNextIteration, but to
					// be explicit here:
					break;
				}
			}
		} catch (IOException e) {
			this.nextByte = -1;
		}

		// At last: finished!
		return this.kmer;
	}

	/**
	 * Returns the reverse complement to the kmer returned by next(). Only
	 * correct if called after next()! Using this is generally cheaper than
	 * calculating the reverse complement by yourself as the reverse complement
	 * is updated byte by byte (i.e. a new byte is read with next(), thus a
	 * single byte of the previous reverse complement needs to be updated) to
	 * avoid duplicate handling of the same byres.
	 *
	 * @return
	 */
	public byte[] getReverseComplement() {
		return this.complement;
	}

	/**
	 * Returns the coordinates of the current k-mer. Only correct if called
	 * after next()!
	 *
	 * @return
	 */
	public KMerCoordinates getCoordinates() {
		return new KMerCoordinates(
				this.recordIndexInFile,
				this.sequenceIndexInFile,
				this.sequenceIndexInRecord,
				this.sequenceIndexInFile + this.skippedKmersInFile,
				this.sequenceIndexInRecord + this.skippedKmersInRecord,
				this.kmer);
	}
}
