package org.husonlab.fmhdist.util;

import jloda.util.FileUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Iterator that produces k-mers of a fixed size for a given FASTA file. This
 * iterator is based on a buffered reader that supplies lines of the file as
 * needed. K-mers will never span multiple sequences in the FASTA.
 */
public class LineKMerIterator implements KMerIterator {
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

	private int k;
	private byte[] kmer;
	private byte[] kmerReverseComplement;
	private byte[] preloadedKmer;
	private byte[] preloadedKmerReverseComplement;

	private byte[] nextLine;
	private byte[] currentLine;
	private int linePointer;

	private boolean isEOF;
	private boolean isPreloaded;
	private BufferedReader reader;

	// Indices to keep track of the origin of the _current_ k-mer. The values
	// will always be copied from the preloaded variants (see below). The values
	// will be fed into the KMerCoordinates if the corresponding method
	// getCoordinates() is called. The separation enables to prepare the next
	// k-mer during a call to "next()" while still be able to load details of
	// the current kmer that is returned by the same call to next().
	private int recordIndexInFile = 0;
	private int skippedKmersInRecord = 0;
	private int sequenceIndexInRecord = 0;
	// Assumption: we call "getCoordinates()" usually less often than "next()",
	// so we save some time
	private int sequencesBeforeRecord = 0;
	private int skippedBeforeRecord = 0;

	// Indices to keep track of the origin of the _next_ k-mer. Those values
	// will actually be incremented as the stream is processed.Those will be
	// written to the current indices during the next call to next().
	private int preloadedRecordIndexInFile = -1;
	private int preloadedSkippedKmersInRecord = 0;
	private int preloadedSequenceIndexInRecord = 0;

	private int preloadedSequencesBeforeRecord = 0;
	private int preloadedSkippedBeforeRecord = 0;

	/**
	 * Creates a new LineKMerIterator that decomposes underlying file into its
	 * k-mers.
	 *
	 * @param k      The k-mer size to apply
	 * @param reader The reader from which the lines of the FASTA file should be
	 *               read
	 * @param skipN  Bool flag, if set to false, k-mers containing ambiguous
	 *               bases won't be skipped.
	 * @throws IOException
	 */
	public LineKMerIterator(int k, BufferedReader reader, boolean skipN) throws IOException {
		this.k = k;
		this.kmer = new byte[k];
		this.kmerReverseComplement = new byte[k];
		this.preloadedKmer = new byte[k];
		this.preloadedKmerReverseComplement = new byte[k];

		this.isEOF = false;
		this.isPreloaded = false;

		this.isAmbiguousChar['N'] = skipN;
		this.isAmbiguousChar['n'] = skipN;

		this.linePointer = 0;
		this.reader = reader;
		String current = reader.readLine();
		String next = reader.readLine();

		if (current == null || next == null) {
			throw new IOException("file is too short, valid FASTA files have at least two lines");
		}
		this.currentLine = current.getBytes();
		this.nextLine = next.getBytes();

		this.preload();
	}

	/**
	 * Creates a new LineKMerIterator that decomposes underlying file into its
	 * k-mers.
	 *
	 * @param k        The k-mer size to apply
	 * @param fileName Path to the Fasta file to read. This could also be a URL
	 *                 or a path to a zip/gzip file.
	 * @param skipN    Bool flag, if set to false, k-mers containing ambiguous
	 *                 bases won't be skipped.
	 * @throws IOException
	 */
	public LineKMerIterator(int k, String fileName, boolean skipN) throws IOException {
		this(k, new BufferedReader(new InputStreamReader(FileUtils.getInputStreamPossiblyZIPorGZIP(fileName))), skipN);
	}

	private void handleSequenceStart() {
		this.preloadedSkippedBeforeRecord += this.preloadedSkippedKmersInRecord;
		this.preloadedSequencesBeforeRecord += this.preloadedSequenceIndexInRecord;

		this.preloadedRecordIndexInFile++;
		this.preloadedSequenceIndexInRecord = 0;
		this.preloadedSkippedKmersInRecord = 0;
	}

	/**
	 * This updates all indices such that they equal the preloaded ones.
	 */
	private void copyIndices() {
		this.sequenceIndexInRecord = this.preloadedSequenceIndexInRecord;
		this.skippedKmersInRecord = this.preloadedSkippedKmersInRecord;
		// TODO: maybe only do when recordIndexInFile changes?
		if (this.recordIndexInFile != this.preloadedRecordIndexInFile) {
			this.skippedBeforeRecord = this.preloadedSkippedBeforeRecord;
			this.sequencesBeforeRecord = this.preloadedSequencesBeforeRecord;
			this.recordIndexInFile = this.preloadedRecordIndexInFile;
		}
	}

	private boolean isSequenceChar() {
		return !isLineContainingSkippableChar[this.currentLine[linePointer]];
	}

	private void feedLine() throws IOException {
		if (this.isEOF) {
			return;
		}

		this.linePointer = 0;
		this.currentLine = this.nextLine;
		String next = reader.readLine();
		if (next == null) {
			this.isEOF = true;
			return;
		}
		this.nextLine = next.getBytes();

		// Peek at the first byte - do we start a new sequence?
		// Assumption: New sequences only occur at the beginning of a new line
		// We can remove the checks for all calls to next()!
		// Only do this if we are not currently preloading
		if (this.currentLine[0] == '>' && !this.isPreloaded) {
			this.currentLine = this.nextLine;
			next = reader.readLine();
			if (next == null) {
				throw new IOException("fasta file contains header without body");
			}
			this.nextLine = next.getBytes();
			if (this.currentLine[0] == '>') {
				throw new IOException("fasta file contains header without body");
			}

			this.handleSequenceStart();
			this.preload();
		}
	}

	private void moveCursor() throws IOException {
		if (++this.linePointer >= this.currentLine.length) {
			this.feedLine();
		}
	}

	private void preload() throws IOException {
		this.isPreloaded = true;
		int i = 0;
		while (this.hasNext()) {
			while (i < this.k - 1 && this.hasNext()) {
				if (isSequenceChar()) {
					if (isAmbiguousChar[this.currentLine[linePointer]]) {
						this.preloadedSkippedKmersInRecord += i + 1;
						i = 0;
						this.moveCursor();
						continue;
					}
					this.preloadedKmer[i] = toUpperTable[this.currentLine[linePointer]];
					this.preloadedKmerReverseComplement[this.k - i - 1] = toUpperTable[complementTable[this.currentLine[linePointer]]];
					this.moveCursor();
					i++;
				} else {
					// no need to check header start explicitely - if NOT a valid sequence character, will skip the current line.
					this.handleSequenceStart();
					this.feedLine();
					i = 0;
				}
			}
			// Now, let's check if the next byte is correct
			if (this.hasNext()) {
				if (this.isSequenceChar()) {
					if (!isAmbiguousChar[this.currentLine[linePointer]]) {
						return;
					} else {
						this.preloadedSkippedKmersInRecord += i + 1;
						this.moveCursor();
						i = 0;
					}
				} else {
					this.handleSequenceStart();
					this.feedLine();
					i = 0;
				}
			}
		}
	}

	@Override
	public boolean hasNext() {
		return !this.isEOF || this.linePointer < this.currentLine.length;
	}

	@Override
	public byte[] next() {
		if (this.isPreloaded) {
			System.arraycopy(this.preloadedKmer, 0, this.kmer, 0, this.k - 1);
			System.arraycopy(this.preloadedKmerReverseComplement, 1, this.kmerReverseComplement, 1, this.k - 1);
			this.isPreloaded = false;
		} else {
			System.arraycopy(this.kmer, 1, this.kmer, 0, this.k - 1);
			System.arraycopy(this.kmerReverseComplement, 0, this.kmerReverseComplement, 1, this.k - 1);
		}
		this.kmer[this.k - 1] = toUpperTable[this.currentLine[linePointer]];
		this.kmerReverseComplement[0] = toUpperTable[complementTable[this.currentLine[linePointer]]];
		this.copyIndices();
		try {
			this.preloadedSequenceIndexInRecord++;
			this.moveCursor();
			if (this.hasNext()) {
				// We only need to check if the next character is ambiguous: if
				// we have a line break and start a new sequence, this is
				// already handled implicitely in the this.moveCursor() -->
				// this.feedLine() chain.
				if (isAmbiguousChar[this.currentLine[this.linePointer]]) {
					this.preloadedSkippedKmersInRecord += this.k;
					this.moveCursor();
					this.isPreloaded = true; // avoid preload in feedLine
					this.preload();
				}
			}
		} catch (IOException e) {
			this.isEOF = true;
			this.linePointer = this.currentLine.length;
		}
		return this.kmer;
	}

	@Override
	public void close() throws IOException {
		this.reader.close();
	}

	@Override
	public int getK() {
		return this.k;
	}

	@Override
	public byte[] getReverseComplement() {
		return this.kmerReverseComplement;
	}

	@Override
	public KMerCoordinates getCoordinates() {
		return new KMerCoordinates(
				this.recordIndexInFile,
				this.sequenceIndexInRecord + this.sequencesBeforeRecord,
				this.sequenceIndexInRecord,
				this.sequenceIndexInRecord + this.sequencesBeforeRecord + this.skippedKmersInRecord + this.skippedBeforeRecord,
				this.sequenceIndexInRecord + this.skippedKmersInRecord,
				this.kmer
		);
	}

}
