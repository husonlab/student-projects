package org.husonlab.fmhdist.util;

import java.io.Closeable;
import java.util.Iterator;

/**
 * A Iterator to support K-Mer decomposition. K-mers returned by this iterator
 * won't span multiple sequences inside the same FASTA file.
 */
public interface KMerIterator extends Iterator<byte[]>, Closeable {
	/**
	 * Returns the k-mer size
	 *
	 * @return
	 */
	public int getK();

	/**
	 * Returns the reverse complement of the k-mer. Callers need to copy the
	 * value as the underlying byte array will be overwritten with the next call
	 * to the iterators "next()" method. The content of the array will always
	 * belong to the last k-mer returned by "next()".
	 *
	 * @return
	 */
	public byte[] getReverseComplement();

	/**
	 * Returns the KMerCoordinates that belong to the last k-mer returned by
	 * "next()"
	 *
	 * @return
	 */
	public KMerCoordinates getCoordinates();
}
