# Contributing to eMERLIN CASA Pipeline

Thank you for considering contributing to the eMERLIN CASA Pipeline! This document provides guidelines and instructions for contributing to this project.

## Code of Conduct

By participating in this project, you agree to abide by our code of conduct. Please be respectful, considerate, and collaborative.

## How Can I Contribute?

### Reporting Bugs

- Check if the bug has already been reported in the [Issues](https://github.com/e-merlin/eMERLIN_CASA_pipeline/issues)
- If not, create a new issue with a clear title and description
- Include steps to reproduce the problem and the expected behavior
- Add relevant information about your environment (OS, Python version, etc.)

### Suggesting Enhancements

- Check if the enhancement has already been suggested in the [Issues](https://github.com/e-merlin/eMERLIN_CASA_pipeline/issues)
- If not, create a new issue with a clear title and description of the proposed enhancement
- Explain why this enhancement would be useful to most users

### Pull Requests

1. Fork the repository and create your branch from `main`
2. If you've added code that should be tested, add tests
3. Ensure your code follows the project's style guidelines
4. Ensure all tests pass
5. Make sure your code lints
6. Issue a pull request

## Development Guidelines

### Environment Setup

We recommend using a conda environment for development:

```bash
git clone https://github.com/e-merlin/eMERLIN_CASA_pipeline.git
cd eMERLIN_CASA_pipeline
conda env create -f environment.yml
conda activate emcp
```

### Project Structure

- `src/eMCP/` - Core source code
- `docs/` - Documentation files
- `tests/` - Test suite

### Code Style

- Follow PEP 8 style guidelines
- Use docstrings for functions and classes
- Keep functions small and focused on a single task
- Use clear, descriptive variable names

### Commit Messages

- Use clear, descriptive commit messages
- Start with a short summary (50 chars or less)
- Reference issues and pull requests where appropriate

### Documentation

- Update documentation when changing code
- Write clear, concise explanations
- Include examples when possible

## Release Process

1. Update version number in `docs/package_version.json`
2. Update CHANGELOG.md with changes since the last release
3. Create a release on GitHub
4. The GitHub Actions workflow will automatically publish to PyPI

Thank you for your contributions!
