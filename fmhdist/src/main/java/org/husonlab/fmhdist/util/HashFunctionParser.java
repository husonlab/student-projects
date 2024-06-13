package org.husonlab.fmhdist.util;

import java.util.ArrayList;
import java.util.List;

import net.openhft.hashing.LongHashFunction;


public class HashFunctionParser {
	public static final String FARM_HASH_NAME = "farm";
	public static final String MURMUR3_HASH_NAME = "murmur3";
	public static final String METRO_HASH_NAME = "metro";
	public static final String XX64_HASH_NAME = "xx64";

	public static List<String> getSupportedFunctions() {
		List<String> result = new ArrayList<>();
		result.add(FARM_HASH_NAME);
		result.add(MURMUR3_HASH_NAME);
		result.add(METRO_HASH_NAME);
		result.add(XX64_HASH_NAME);
		return result;
	}

	public static LongHashFunction createHashFunction(String name, long seed) throws IllegalArgumentException {
		switch (name) {
			case FARM_HASH_NAME:
				return LongHashFunction.farmUo(seed);
			case MURMUR3_HASH_NAME:
				return LongHashFunction.murmur_3(seed);
			case METRO_HASH_NAME:
				return LongHashFunction.metro(seed);
			case XX64_HASH_NAME:
				return LongHashFunction.xx(seed);
		}

		throw new IllegalArgumentException("unkown hash function name");
	}
}
