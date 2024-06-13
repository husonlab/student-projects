package org.husonlab.fmhdist.util.experimental;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import org.husonlab.fmhdist.util.KMerCoordinates;
import org.husonlab.fmhdist.util.KMerIterator;

import jloda.util.FileLineIterator;

public class InMemoryKMerIterator implements KMerIterator {
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

	private List<byte[]> sequences;
	private int sequencePointer = 0;

	private byte[] kmer;
	private byte[] complement;

	private byte[] currentSequence;

	private int charPointer = 0;

	private int groupHelper;

	public InMemoryKMerIterator(int k, String fileName, boolean skipN) throws IOException {
		this.k = k;

		this.kmer = new byte[k];
		this.complement = new byte[k];

		this.isAmbiguousChar['N'] = skipN;
		this.isAmbiguousChar['n'] = skipN;

		synchronized (InMemoryKMerIterator.class) {
			try (FileLineIterator it = new FileLineIterator(fileName)) {
				this.sequences = it
						.stream()
						.collect(Collectors.groupingBy(line -> {
							if (line.startsWith(">")) {
								return ++groupHelper;
							}
							return groupHelper;
						})).entrySet()
						.stream()
						.map(l -> l.getValue()
								.stream()
								.filter(line -> !line.startsWith(">"))
								.map(line -> line.replaceAll("\\s+", "").toUpperCase())
								.collect(Collectors.joining())
								.getBytes())
						.collect(Collectors.toList());
			}
		}

		this.currentSequence = sequences.get(this.sequencePointer++);
		this.preload();
	}

	private void feedSequence() {
		if (this.sequencePointer >= this.sequences.size()) {
			return;
		}

		this.charPointer = 0;
		this.currentSequence = this.sequences.get(this.sequencePointer++);
	}

	/**
	 * Return true if new sequence was loaded.
	 *
	 * @return
	 * @throws IOException
	 */
	private boolean moveCursor() {
		if (++this.charPointer >= this.currentSequence.length) {
			this.feedSequence();
			return true;
		}
		return false;
	}

	private void preload() {
		while (true) {
			int i = 0;
			while (i < this.k - 1 && this.hasNext()) {
				// No need to check if is sequence char - we are only handling
				// sequence characters!
				if (isAmbiguousChar[this.currentSequence[charPointer]]) {
					this.moveCursor();
					i = 0;
					continue;
				}
				if (this.moveCursor()) {
					i = 0;
					continue;
				}
				i++;
			}

			// Peek at the next byte, is it valid?
			// Do not move char pointer if it is valid!
			if (this.hasNext()) {
				if (!isAmbiguousChar[this.currentSequence[this.charPointer]]) {
					return;
				}
			} else {
				return;
			}
		}
	}

	@Override
	public boolean hasNext() {
		// charPoint <= currentSequence.length as this always points at the
		// _next_ byte (i.e. the k+1 byte for the last k-mer)
		return this.sequencePointer < this.sequences.size() || this.charPointer < this.currentSequence.length;
	}

	@Override
	public byte[] next() {
		System.arraycopy(this.currentSequence, this.charPointer - this.k + 1, this.kmer, 0, this.k);
		for (int i = 0; i < this.k; i++) {
			this.complement[i] = complementTable[this.kmer[this.k - i - 1]];
		}

		if (this.moveCursor() || isAmbiguousChar[this.currentSequence[this.charPointer]]) {
			this.preload();
		}

		return this.kmer;
	}

	@Override
	public void close() throws IOException {
		return;
	}

	@Override
	public int getK() {
		return this.k;
	}

	@Override
	public byte[] getReverseComplement() {
		return this.complement;
	}

	/**
	 * Won't return actual coordinates, maybe I can split the interface?
	 */
	@Override
	public KMerCoordinates getCoordinates() {
		return new KMerCoordinates(0, 0, 0, 0, 0, this.kmer);
	}

}
