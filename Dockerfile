FROM astral/uv:python3.12-bookworm-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    build-essential \
    graphviz \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Copy project files
COPY pyproject.toml uv.lock* ./
COPY src/ ./src/

# Install dependencies
RUN uv sync --no-dev

# Add venv to PATH so CLI commands are available
ENV PATH="/app/.venv/bin:$PATH"

# Set entrypoint
CMD ["bash"]
