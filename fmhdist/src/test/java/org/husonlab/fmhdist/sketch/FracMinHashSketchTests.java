package org.husonlab.fmhdist.sketch;

import net.openhft.hashing.LongHashFunction;
import org.husonlab.fmhdist.util.KMerCoordinates;
import org.husonlab.fmhdist.util.KMerIterator;
import org.husonlab.fmhdist.util.experimental.FastKMerIterator;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

public class FracMinHashSketchTests {
	class TestKMerIterator implements KMerIterator {
		String[] kmers = new String[]{"ATGC", "TGCA", "GCAA", "CAAG"};
		String[] complements = new String[]{"GCAT", "TGCA", "TTGC", "CTTG"};
		int i = -1;

		@Override
		public boolean hasNext() {
			return i < kmers.length - 1;
		}

		@Override
		public byte[] next() {
			return kmers[++i].getBytes();
		}

		@Override
		public void close() throws IOException {
			return;
		}

		@Override
		public int getK() {
			return 4;
		}

		@Override
		public byte[] getReverseComplement() {
			return complements[i].getBytes();
		}

		@Override
		public KMerCoordinates getCoordinates() {
			return new KMerCoordinates(0, 0, 0, 0, 0, new byte[1]);
		}
	}

	/**
	 * @throws IOException
	 */
	@Test
	public void shouldCalculateFracMinSketch() throws IOException {
		try (FastKMerIterator kmers = new FastKMerIterator(21, "src/test/resources/virus1.fasta", true)) {
			FracMinHashSketch sketch = FracMinHashSketch.compute("test", kmers, 10, LongHashFunction.murmur_3(42), 42, true);

			// the test file contains the kmer "GAACTTTGACTGGTTTTTGGG". For
			// this, the reverse complement is "CCCAAAAACCAGTCAAAGTTC", which
			// is smaller.
			long expectedHash = LongHashFunction.murmur_3(42).hashBytes("CCCAAAAACCAGTCAAAGTTC".getBytes());

			// The expected hash is -8556973683090215967, which is lesser than
			// Long.MIN + (Long.MAX *  1/s * 2) with s = 10 as used above.
			// Thus, we expect it to be included in the result.

			List<Long> values = new ArrayList<>();
			for (long value : sketch.getValues()) {
				values.add(value);
			}
			assertThat(values, hasItems(
					expectedHash
			));
			System.out.println();
		}
	}

	@Test
	public void differentHashingFunctionsShouldGenerateIncompatibleSerializations() throws IOException {
		KMerIterator it1 = new TestKMerIterator();
		KMerIterator it2 = new TestKMerIterator();
		byte[] string1 = FracMinHashSketch.compute("Test", it1, 1, LongHashFunction.farmNa(42), 42, false).getBytes();
		byte[] string2 = FracMinHashSketch.compute("Test", it2, 1, LongHashFunction.murmur_3(42), 42, false).getBytes();

		FracMinHashSketch sketch1 = FracMinHashSketch.parse(string1);
		FracMinHashSketch sketch2 = FracMinHashSketch.parse(string2);

		assertThat(sketch1.getHashedMagicNumber(), not(equalTo(sketch2.getHashedMagicNumber())));
	}

	@Test
	public void sametHashingFunctionsShouldGenerateCompatibleSerializations() throws IOException {
		KMerIterator it1 = new TestKMerIterator();
		KMerIterator it2 = new TestKMerIterator();
		byte[] string1 = FracMinHashSketch.compute("Test", it1, 1, LongHashFunction.farmNa(42), 42, false).getBytes();
		byte[] string2 = FracMinHashSketch.compute("Test", it2, 1, LongHashFunction.farmNa(42), 42, false).getBytes();

		FracMinHashSketch sketch1 = FracMinHashSketch.parse(string1);
		FracMinHashSketch sketch2 = FracMinHashSketch.parse(string2);

		assertThat(sketch1.getHashedMagicNumber(), equalTo(sketch2.getHashedMagicNumber()));
	}


	@Test
	@Ignore // Not a real test, but for understanding the performance
	public void benchmarkBloomPerformance() throws IOException {
		String url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/225/685/GCA_026225685.1_Pinf1306_UCR_1/GCA_026225685.1_Pinf1306_UCR_1_genomic.fna.gz";
		int k = 21;
		int s = 1000;
		long initialHeap = Runtime.getRuntime().totalMemory();
		long start = System.currentTimeMillis();

		try (FastKMerIterator kmers = new FastKMerIterator(k, url, true)) {
			FracMinHashSketch.compute(
					"test",
					kmers,
					s,
					LongHashFunction.farmUo(42),
					42,
					true
			);
			long finalHeap = Runtime.getRuntime().totalMemory();
			long end = System.currentTimeMillis();
			System.out.println(String.format("runtime: %d", (end - start) / 1000));
			System.out.println(String.format("initial heap: %d\nfinal heap: %d", initialHeap / 1024 / 1024, finalHeap / 1024 / 1024));
		}
	}

	// @Ignore // Not a real test, but for understanding the performance
	@Test
	public void testSegments() throws IOException {
		String url = "/home/gnampfelix/Documents/term5/thesis/data/phytophthora/test_genomes/PhyInf1.fna.gz";
		long initialHeap = Runtime.getRuntime().totalMemory();
		long start = System.currentTimeMillis();
		try (FastKMerIterator kmers = new FastKMerIterator(21, url, true)) {
			FracMinHashSketch sketch = FracMinHashSketch.compute("test", kmers, 1000, LongHashFunction.farmUo(42), 42, true);
			long finalHeap = Runtime.getRuntime().totalMemory();
			long end = System.currentTimeMillis();
			System.out.println(String.format("sketch size: %d", sketch.getValues().length));
			System.out.println(String.format("runtime: %d", (end - start) / 1000));
			System.out.println(String.format("initial heap: %d\nfinal heap: %d", initialHeap / 1024 / 1024, finalHeap / 1024 / 1024));
		}
	}

}
