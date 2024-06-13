package org.husonlab.fmhdist.ncbi;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;
import org.openapitools.client.ApiException;

public class NcbiApiTests {
	@Test
	public void getGenomeWithLink() throws ApiException, AmbiguousDataException {
		NcbiApi api = new NcbiApi("https://api.ncbi.nlm.nih.gov/datasets/v2alpha");
		List<String> accessionCodes = new ArrayList<>(Arrays.asList("GCF_000819615.1"));
		List<Genome> genomes = api.getGenomes(accessionCodes);
		assertEquals(1, genomes.size());
		assertEquals("Escherichia phage phiX174", genomes.get(0).getOrganismName());
		assertNotNull(genomes.get(0).getDownloadLink());
	}

	@Test
	public void getGenomesWithLink() throws ApiException, AmbiguousDataException {
		NcbiApi api = new NcbiApi("https://api.ncbi.nlm.nih.gov/datasets/v2alpha");
		List<String> accessionCodes = new ArrayList<>(Arrays.asList("GCF_000819615.1", "GCF_000836805.1"));
		List<Genome> genomes = api.getGenomes(accessionCodes);
		assertEquals(2, genomes.size());
		assertNotNull(genomes.get(0).getDownloadLink());
		assertNotNull(genomes.get(1).getDownloadLink());
	}

	@Test
	public void getInvalidAccessionCode() throws ApiException, AmbiguousDataException {
		NcbiApi api = new NcbiApi("https://api.ncbi.nlm.nih.gov/datasets/v2alpha");
		List<String> accessionCodes = new ArrayList<>(Arrays.asList("GF_000819615.1"));
		List<Genome> genomes = api.getGenomes(accessionCodes);
		assertEquals(0, genomes.size());
	}

	@Test
	public void getTaxonomy() throws ApiException, AmbiguousDataException {
		NcbiApi api = new NcbiApi("https://api.ncbi.nlm.nih.gov/datasets/v2alpha");
		List<String> accessionCodes = new ArrayList<>(Arrays.asList("GCF_000819615.1", "GCF_000836805.1"));
		List<Genome> genomes = api.getGenomes(accessionCodes);
		TaxonomyTree tree = api.getTaxonomyTreeForGenomes(genomes);

		// the tree that has those two leaves has 16 nodes in total
		assertEquals(16, tree.getTree().nodes().size());
	}


}
