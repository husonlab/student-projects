package org.husonlab.fmhdist.sketch;

import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class SequenceGrouperTests {
	@Test
	public void shouldGroupSequences() throws IOException {
		List<byte[]> result = new SequenceGrouper("src/test/resources/virus2.fasta").getGroups();
		assertEquals(2, result.size());
		assertEquals("abcdefghijkl", new String(result.get(0)));
		assertEquals("1234567890", new String(result.get(1)));
	}

	@Test
	public void shouldGroupSingleSequences() throws IOException {
		List<byte[]> result = new SequenceGrouper("src/test/resources/virus1.fasta").getGroups();
		assertEquals(1, result.size());
		assertEquals("ACTCGTATGAACTTTGACTGGTTTTTGGGGCGCGAGAGTTTGGTTTGGATGAAACGCACCCGCTATCGAT" +
					 "AGATAGTTTACCTTTCTCTTTTTTTGAAAACGCACCCGCTATCGATAGATAGGGTTACCCTTTCTCTTTG" +
					 "AGGTTTTTGGGAAAACGCACCTGCTATCGATAGTTTAACCTTTCTCTTTGAGGTTTTTGGATGAAACGCA" +
					 "CCCGCTATCTATCGTTACACTTCTCTATCGCATTAACTTGTCTTTTGTAACTTTGCTCGTTTTGTTAGGA" +
					 "AAGATAGATACACGAAGGAAGCACTCGCTCGTTAGGAAAGATACACACGAGGAAGCACTCGCTCGTTAGG" +
					 "AAAGATAGATACACGAGGAAGCCTCCCTCGCTAGACACTAGGAAGTCATCGCTCGTTAGGAAAGACATAA", new String(result.get(0)));
	}
}
