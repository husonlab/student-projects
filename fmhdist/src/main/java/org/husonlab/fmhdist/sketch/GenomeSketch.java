package org.husonlab.fmhdist.sketch;

import java.io.IOException;
import java.util.logging.Logger;

import org.husonlab.fmhdist.ncbi.Genome;
import org.husonlab.fmhdist.util.KMerIterator;
import org.husonlab.fmhdist.util.LineKMerIterator;

import net.openhft.hashing.LongHashFunction;

/**
 * Composite to connect a sketch to a given genome.
 */
public class GenomeSketch {
	private Genome genome;
	private FracMinHashSketch sketch;
	private static Logger logger = Logger.getLogger(GenomeSketch.class.getName());

	/**
	 * Creates a new GenomeSketch by calculating the FracMinHash sketch of the
	 * given genome.
	 *
	 * @param genome             the input genome to sketch
	 * @param kSize              the k-mer size to apply
	 * @param sParam             the scaling parameter to apply
	 * @param hashFunction       the hash function to use
	 * @param seed               the random seed that was used to generate the hash function
	 * @param prepareCoordinates boolean flag to indicate if the KMerCoordinates
	 *                           will be exported later.
	 * @return A new GenomeSketch
	 * @throws IOException
	 */
	public static GenomeSketch sketch(Genome genome, int kSize, int sParam, LongHashFunction hashFunction, int seed, boolean prepareCoordinates) throws IOException {
		logger.fine("Calculating sketch for " + genome.getAccession());
		final GenomeSketch result = new GenomeSketch(genome);

		KMerIterator kmers;
		kmers = new LineKMerIterator(kSize, genome.getFastaUrl(), true);
		result.sketch = FracMinHashSketch.compute(genome.getAccession(), kmers, sParam, hashFunction, seed, prepareCoordinates);
		kmers.close();
		return result;
	}

	private GenomeSketch(Genome genome) {
		this.genome = genome;
	}

	public GenomeSketch(Genome genome, FracMinHashSketch sketch) {
		this.genome = genome;
		this.sketch = sketch;
	}

	public Genome getGenome() {
		return this.genome;
	}

	public FracMinHashSketch getSketch() {
		return this.sketch;
	}
}
