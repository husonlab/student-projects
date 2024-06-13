package org.husonlab.fmhdist.util.experimental;

import java.io.BufferedReader;
import java.io.IOException;

public interface FileConsumer {
	public boolean isReady();

	public void setReady(boolean v);

	public String getLine() throws IOException;

	public void setReader(BufferedReader reader);
}
