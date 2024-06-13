package org.husonlab.fmhdist.sketch;

import net.openhft.hashing.LongHashFunction;
import org.husonlab.fmhdist.util.experimental.FastKMerIterator;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.equalTo;

public class DistanceTests {
	@Test
	public void testIntersectionSize() {
		long[] a = new long[]{1, 2, 3};
		long[] b = new long[]{2};
		assertThat(Distance.getIntersectionSize(a, b), equalTo(1));
	}

	@Test
	public void testJaccardIndex() {
		long[] a = new long[]{1, 2, 3};
		long[] b = new long[]{2};
		assertThat(Distance.calculateJaccardIndex(a, b, 1), equalTo(1.0 / 3.0));
	}

	@Test
	public void testJaccardIndexWithSelf() {
		long[] a = new long[]{1, 2, 3};
		assertThat(Distance.calculateJaccardIndex(a, a, 1), equalTo(1.0));
	}

	@Test
	public void testContainmentIndex() {
		long[] a = new long[]{1, 2, 3};
		long[] b = new long[]{2};
		assertThat(Distance.calculateContainmentIndex(b, a, 1), equalTo(1.0));
	}

	@Test
	public void testJaccrdDistanceWithSelf() {
		assertThat(Distance.jaccardToDistance(1, 21), equalTo(0.0));
	}

	// Not really a test, more an example!
	@Test
	@Ignore
	public void runExampleDistanceCalculation() throws IOException, IncompatibleParameterException {
		final String[] downloadUrls = new String[]{
				"/home/gnampfelix/Documents/term5/thesis/data/phytophthora/test_genomes/EColi1.fna.gz",
				"/home/gnampfelix/Documents/term5/thesis/data/phytophthora/test_genomes/EColi2.fna.gz",
				"/home/gnampfelix/Documents/term5/thesis/data/phytophthora/test_genomes/PhyAga.fna.gz",
				"/home/gnampfelix/Documents/term5/thesis/data/phytophthora/test_genomes/PhyInf1.fna.gz",
				"/home/gnampfelix/Documents/term5/thesis/data/phytophthora/test_genomes/PhyInf2.fna.gz",
				"/home/gnampfelix/Documents/term5/thesis/data/phytophthora/test_genomes/PhyKer.fna.gz"
		};

		final String[] taxa = new String[]{
				"E. coli",
				"E. coli",
				"P. agathidic",
				"P. infestans",
				"P. infestans",
				"P. kernoviae"
		};

		long start = System.currentTimeMillis();
		final int k = 21;
		final int s = 1000;
		List<FracMinHashSketch> sketches = new ArrayList<>();
		for (String url : downloadUrls) {
			try (FastKMerIterator kmers = new FastKMerIterator(k, url, true)) {
				sketches.add(FracMinHashSketch.compute(
						url,
						kmers,
						s,
						LongHashFunction.farmUo(42),
						42,
						true
				));
			}
			System.out.println(String.format("%d seconds", (System.currentTimeMillis() - start) / 1000));
		}

		// each entry consists of the three different distances
		List<List<List<Double>>> distances = new ArrayList<>();
		for (int i = 0; i < sketches.size(); i++) {
			List<List<Double>> currentRow = new ArrayList<>();
			for (int j = i; j < sketches.size(); j++) {
				List<Double> currentDistances = new ArrayList<>();
				double jaccard = Distance.calculateJaccardIndex(sketches.get(i).getValues(),
						sketches.get(j).getValues(), s);
				currentDistances.add(Distance.jaccardToDistance(jaccard, k));
				currentDistances.add(jaccard);
				currentRow.add(currentDistances);
			}
			distances.add(currentRow);
		}

		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < distances.size(); j++) {
				System.out.print(String.format("%s:\t", taxa[j]));
				for (int l = 0; l < distances.get(j).size(); l++) {
					System.out.print(String.format("%f\t", distances.get(j).get(l).get(i)));
				}
				System.out.print("\n");
			}
			System.out.print("\n");
		}
		System.out.println();

	}
}
