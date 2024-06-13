package org.husonlab.fmhdist.util;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.equalTo;

import org.junit.Test;

public class KMerCoordinatesTest {
	@Test
	public void itShouldParseCoordinates() {
		String input = "2901,203233312,1374,203233312,1374,ACCATTACAATGACTTTGGAT,17,0";
		KMerCoordinates coords = KMerCoordinates.fromString(input);
		assertThat(coords.getRecordIndexInFile(), equalTo(2901));
		assertThat(coords.getSequenceIndexInFile(), equalTo(203233312));
		assertThat(coords.getSequenceIndexInRecord(), equalTo(1374));
		assertThat(coords.getSequenceIndexInFileIncludingAmbiguous(), equalTo(203233312));
		assertThat(coords.getSequenceIndexInRecordIncludingAmbiguous(), equalTo((1374)));
		assertThat(coords.getKmer(), equalTo("ACCATTACAATGACTTTGGAT".getBytes()));
	}
}
