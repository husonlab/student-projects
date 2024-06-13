package org.husonlab.fmhdist.util.experimental;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import org.husonlab.fmhdist.util.KMerCoordinates;
import org.husonlab.fmhdist.util.KMerIterator;

import jloda.util.FileUtils;

public class LineKMerIteratorConsumer implements KMerIterator {
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

	private FileConsumer consumer;

	public LineKMerIteratorConsumer(int k, String fileName, FileProducer producer, int threadIndex, boolean skipN)
			throws IOException {
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

		this.reader = new BufferedReader(new InputStreamReader(FileUtils.getInputStreamPossiblyZIPorGZIP(fileName)));
		this.consumer = producer.getFileConsumer(threadIndex);
		this.consumer.setReader(this.reader);
		this.consumer.setReady(false);

		String current;
		String next;

		// the consumer should be ready semantically because it is set to ready
		// before above code can obtain it. Still, let's ensure this here.
		while (true) {
			if (this.consumer.isReady()) {
				current = this.consumer.getLine();
				this.consumer.setReady(false);
				break;
			}
		}

		while (true) {
			if (this.consumer.isReady()) {
				next = this.consumer.getLine();
				this.consumer.setReady(false);
				break;
			}
		}

		if (current == null || next == null) {
			throw new IOException("file is too short, valid FASTA files have at least two lines");
		}
		this.currentLine = current.getBytes();
		this.nextLine = next.getBytes();

		this.preload();
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

		String next;
		while (true) {
			if (this.consumer.isReady()) {
				next = this.consumer.getLine();
				if (next == null) {
					this.isEOF = true;
					// don't un-ready the reader now - we want to stay in control
					return;
				}
				break;
			}
		}

		this.consumer.setReady(false);
		this.nextLine = next.getBytes();

		// Peek at the first byte - do we start a new sequence?
		// Assumption: New sequences only occur at the beginning of a new line
		// We can remove the checks for all calls to next()!
		// Only do this if we are not currently preloading
		if (this.currentLine[0] == '>' && !this.isPreloaded) {
			this.currentLine = this.nextLine;
			while (true) {
				if (this.consumer.isReady()) {
					next = this.consumer.getLine();
					// only ready the consumer if we are sure that we want to read afterwards
					if (next == null) {
						throw new IOException("fasta file contains header without body");
					}
					if (this.currentLine[0] == '>') {
						throw new IOException("fasta file contains header without body");
					}
					this.consumer.setReady(false);
					break;
				}
			}
			this.nextLine = next.getBytes();

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
						i = 0;
						this.moveCursor();
						continue;
					}
					this.preloadedKmer[i] = toUpperTable[this.currentLine[linePointer]];
					this.preloadedKmerReverseComplement[this.k - i
														- 1] = toUpperTable[complementTable[this.currentLine[linePointer]]];
					this.moveCursor();
					i++;
				} else {
					// no need to check header start explicitely - if NOT a valid sequence
					// character, will skip the current line.
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
						this.moveCursor();
						i = 0;
					}
				} else {
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
		try {
			this.moveCursor();
			if (this.hasNext()) {
				// We only need to check if the next character is ambiguous: if
				// we have a line break and start a new sequence, this is
				// already handled implicitely in the this.moveCursor() -->
				// this.feedLine() chain.
				if (isAmbiguousChar[this.currentLine[this.linePointer]]) {
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
		while (true) {
			if (this.consumer.isReady()) {
				// ensure that the producer is not trying to access the reader
				this.reader.close();
				return;
			}
		}
	}

	@Override
	public int getK() {
		return this.k;
	}

	@Override
	public byte[] getReverseComplement() {
		return this.kmerReverseComplement;
	}

	/**
	 * Won't return actual coordinates, maybe I can split the interface?
	 */
	@Override
	public KMerCoordinates getCoordinates() {
		return new KMerCoordinates(0, 0, 0, 0, 0, this.kmer);
	}

}
