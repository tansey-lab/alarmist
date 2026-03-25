.PHONY: dev install-uv install-java install-nextflow install-nf-test install-deps install-hooks test test-nf lint format

# Development environment setup
dev: install-uv install-java install-nextflow install-nf-test install-deps install-hooks
	@echo "Development environment setup complete!"

# Run pytest
test:
	uv run pytest

# Run Nextflow integration tests
test-nf:
	@echo "Running Nextflow integration tests..."
	cd nextflow && nextflow run main.nf -profile test

# Run linting
lint:
	uv run ruff check .

# Run formatting
format:
	uv run ruff format .

install-uv:
	@if command -v uv >/dev/null 2>&1; then \
		echo "uv is already installed"; \
	else \
		echo "Installing uv..."; \
		curl -LsSf https://astral.sh/uv/install.sh | sh; \
	fi

install-java:
	@echo "Checking for Java..."
	@if which java >/dev/null 2>&1; then \
		echo "Java is already installed: $$(java -version 2>&1 | head -1)"; \
	else \
		echo "Java not found, installing Temurin JRE 21 from GitHub..."; \
		mkdir -p ~/.local/java; \
		OS=$$(uname -s | tr '[:upper:]' '[:lower:]'); \
		ARCH=$$(uname -m); \
		if [ "$$ARCH" = "x86_64" ]; then ARCH="x64"; fi; \
		if [ "$$ARCH" = "arm64" ]; then ARCH="aarch64"; fi; \
		if [ "$$OS" = "darwin" ]; then OS="mac"; fi; \
		echo "Detected OS=$$OS, ARCH=$$ARCH"; \
		DOWNLOAD_URL=$$(curl -sL https://api.github.com/repos/adoptium/temurin21-binaries/releases/latest \
			| grep "browser_download_url.*jre_$${ARCH}_$${OS}_.*tar.gz\"" \
			| grep -v ".sig" \
			| grep -v "alpine" \
			| head -1 \
			| cut -d'"' -f4); \
		echo "Downloading from: $$DOWNLOAD_URL"; \
		if [ -z "$$DOWNLOAD_URL" ]; then echo "ERROR: Could not find JRE download URL for OS=$$OS ARCH=$$ARCH"; exit 1; fi; \
		curl -fSL "$$DOWNLOAD_URL" -o /tmp/jre21.tar.gz; \
		tar xzf /tmp/jre21.tar.gz -C ~/.local/java --strip-components=1; \
		rm /tmp/jre21.tar.gz; \
		ln -sf ~/.local/java/bin/java ~/.local/bin/java; \
		echo "Java 21 JRE installed to ~/.local/java (symlinked to ~/.local/bin/java)"; \
	fi

install-nextflow:
	@if command -v nextflow >/dev/null 2>&1; then \
		echo "Nextflow is already installed"; \
	else \
		echo "Installing Nextflow..."; \
		mkdir -p ~/.local/bin; \
		NF_VERSION=$$(curl -sL https://api.github.com/repos/nextflow-io/nextflow/releases/latest | grep '"tag_name"' | cut -d'"' -f4); \
		curl -sL "https://github.com/nextflow-io/nextflow/releases/download/$${NF_VERSION}/nextflow" -o ~/.local/bin/nextflow; \
		chmod +x ~/.local/bin/nextflow; \
		echo "Nextflow $${NF_VERSION} installed to ~/.local/bin"; \
	fi

install-nf-test:
	@if command -v nf-test >/dev/null 2>&1; then \
		echo "nf-test is already installed"; \
	else \
		echo "Installing nf-test..."; \
		mkdir -p ~/.local/bin ~/.nf-test && \
		NFT_VERSION=$$(curl -sL https://api.github.com/repos/askimed/nf-test/releases/latest | grep '"tag_name"' | cut -d'"' -f4 | sed 's/^v//') && \
		curl -sL "https://github.com/askimed/nf-test/releases/download/v$${NFT_VERSION}/nf-test-$${NFT_VERSION}.tar.gz" | tar xz -C ~/.nf-test && \
		ln -sf ~/.nf-test/nf-test ~/.local/bin/nf-test && \
		JAVA_CMD=$$(which java) ~/.local/bin/nf-test version && \
		echo "nf-test $${NFT_VERSION} installed to ~/.local/bin"; \
	fi

install-deps:
	@echo "Installing Python dependencies with uv..."
	uv sync --group dev

install-hooks:
	@echo "Installing pre-commit hooks..."
	uv run pre-commit install
	@echo "Pre-commit hooks installed and active!"
