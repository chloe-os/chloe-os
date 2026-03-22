# Contributing

We welcome contributions from bioinformaticians, veterinary researchers, software engineers, and anyone passionate about pet health.

## Getting Started

```bash
# Clone the repo
git clone https://github.com/chloe-os/chloe-os.git
cd chloe-os

# Install dependencies (including dev tools)
uv sync --all-extras

# Run tests
uv run pytest

# Run linter
uv run ruff check .

# Run type checker
uv run mypy packages/chloe-core/src packages/chloe-cli/src
```

## How to Contribute

### Report Issues

Found a bug or have a feature request? [Open an issue](https://github.com/chloe-os/chloe-os/issues) on GitHub.

### Submit Pull Requests

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Make your changes
4. Add tests for new functionality
5. Run the test suite (`uv run pytest`)
6. Run the linter (`uv run ruff check .`)
7. Submit a pull request

### Areas Where We Need Help

- **Canine DLA allele databases** — Better breed-specific DLA allele frequencies
- **Validation studies** — Comparing predictions against known immunogenic neoantigens
- **Additional species support** — Feline, equine, etc.
- **Documentation** — Improving guides for pet owners and researchers
- **Testing** — More comprehensive test cases with real-world data
- **MHC predictor backends** — NetMHCpan integration, custom canine models

## Code Style

- Python 3.11+ with type hints
- Formatted with `ruff format`
- Linted with `ruff check`
- Type-checked with `mypy --strict`

## Project Principles

1. **Accessibility first** — The tool should be usable by pet owners with no bioinformatics background
2. **Transparency** — All algorithms and scoring methods must be documented and open
3. **Safety** — Prominent disclaimers at every stage; never suggest replacing veterinary advice
4. **Modularity** — Each pipeline stage should be independently usable and testable
