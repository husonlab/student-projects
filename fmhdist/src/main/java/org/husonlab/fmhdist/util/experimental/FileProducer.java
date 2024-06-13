package org.husonlab.fmhdist.util.experimental;

import java.io.BufferedReader;
import java.io.IOException;

public class FileProducer {
	private class Consumer implements FileConsumer {
		private volatile String currentLine;
		private volatile String nextLine;
		private volatile BufferedReader reader;
		private volatile boolean isReady;
		private volatile IOException e;
		private volatile boolean hasError;
		private boolean isInitialized;

		public Consumer() {

		}

		public void setReader(BufferedReader reader) {
			this.reader = reader;
			this.isInitialized = false;
			this.hasError = false;
			this.currentLine = null;
			this.nextLine = null;
		}

		/**
		 * true if ready to be used by consumer. Set to false if ready to be
		 * used by producer.
		 */
		@Override
		public boolean isReady() {
			return this.isReady;
		}

		@Override
		public void setReady(boolean v) {
			this.isReady = v;
		}

		public void feedLine() {
			this.currentLine = this.nextLine;
			try {
				this.nextLine = reader.readLine();
			} catch (IOException e) {
				this.hasError = true;
				this.e = e;
			}
		}

		@Override
		public String getLine() throws IOException {
			if (this.hasError) {
				throw this.e;
			}
			return this.currentLine;
		}
	}

	// Caution: need to check thread safety
	private Consumer[] consumers;
	private volatile boolean isClosed;

	public FileProducer(int nThreads) {
		this.consumers = new Consumer[nThreads];
		for (int i = 0; i < nThreads; i++) {
			this.consumers[i] = new Consumer();
			// Consumers are ready initially such that we don't accidentally try
			// to read from them.
			this.consumers[i].isReady = true;
		}
	}

	public FileConsumer getFileConsumer(int threadIndex) throws IOException {
		return this.consumers[threadIndex];
	}

	public void run() {
		while (true) {
			for (int i = 0; i < this.consumers.length; i++) {
				if (!this.consumers[i].isReady) {
					// this thread now has exclusive access
					this.consumers[i].feedLine();
					if (!this.consumers[i].isInitialized) {
						this.consumers[i].feedLine();
						this.consumers[i].isInitialized = true;
					}
					this.consumers[i].isReady = true;
				}
			}
			if (this.isClosed) {
				return;
			}
		}
	}

	public void close() {
		this.isClosed = true;
	}
}
