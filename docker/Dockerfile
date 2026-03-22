FROM python:3.11-slim

LABEL maintainer="Chloe OS Contributors"
LABEL description="Chloe OS — Personalized cancer vaccine pipeline for pets"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    perl \
    cpanminus \
    build-essential \
    libdbi-perl \
    libdbd-mysql-perl \
    libwww-perl \
    libjson-perl \
    curl \
    wget \
    unzip \
    tabix \
    && rm -rf /var/lib/apt/lists/*

# Install Ensembl VEP
ENV VEP_VERSION=113
RUN cd /opt && \
    git clone --depth 1 --branch release/${VEP_VERSION} \
    https://github.com/Ensembl/ensembl-vep.git && \
    cd ensembl-vep && \
    perl INSTALL.pl --AUTO a --NO_HTSLIB --NO_UPDATE 2>&1 | tail -5

ENV PATH="/opt/ensembl-vep:${PATH}"

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# Set up working directory
WORKDIR /app

# Copy project files
COPY pyproject.toml .
COPY packages/ packages/

# Install Python dependencies
RUN uv sync --no-dev

# Install MHCflurry and download models
RUN uv pip install mhcflurry && \
    uv run mhcflurry-downloads fetch models_class1_presentation

# Copy remaining files
COPY . .

# Default command
ENTRYPOINT ["uv", "run", "chloe"]
CMD ["--help"]
