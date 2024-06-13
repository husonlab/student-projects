package org.husonlab.fmhdist.sketch;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import jloda.util.FileLineIterator;

public class SequenceGrouper {
	private int group;
	private List<byte[]> result;

	public SequenceGrouper(String path) throws IOException {
		try (FileLineIterator it = new FileLineIterator(path)) {
			this.result = it
					.stream()
					.collect(Collectors.groupingBy(line -> {
						if (line.startsWith(">")) {
							return ++group;
						}
						return group;
					})).entrySet()
					.stream()
					.map(l -> l.getValue()
							.stream()
							.filter(line -> !line.startsWith(">"))
							.map(line -> line.replaceAll("\\s+", ""))
							.collect(Collectors.joining())
							.getBytes())
					.collect(Collectors.toList());
		}
	}

	public List<byte[]> getGroups() {
		return this.result;
	}
}
