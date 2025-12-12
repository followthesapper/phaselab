# Contributing to PhaseLab

Thank you for your interest in contributing to PhaseLab! This document provides guidelines and instructions for contributing.

## Getting Started

### Prerequisites

- Python 3.9 or higher
- Git
- A GitHub account

### Development Setup

1. Fork the repository on GitHub

2. Clone your fork:
```bash
git clone https://github.com/YOUR_USERNAME/phaselab.git
cd phaselab
```

3. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # Linux/macOS
# or
venv\Scripts\activate  # Windows
```

4. Install development dependencies:
```bash
pip install -e ".[dev]"
```

5. Verify the setup:
```bash
pytest
```

## Development Workflow

### Creating a Branch

Create a branch for your work:

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/your-bug-fix
```

Branch naming conventions:
- `feature/` - New features
- `fix/` - Bug fixes
- `docs/` - Documentation changes
- `refactor/` - Code refactoring
- `test/` - Test additions or fixes

### Making Changes

1. Write your code following the style guidelines below
2. Add tests for new functionality
3. Update documentation if needed
4. Run tests locally before committing

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=phaselab --cov-report=term-missing

# Run specific test file
pytest tests/test_crispr.py

# Run specific test
pytest tests/test_crispr.py::test_guide_design
```

### Code Style

We use the following tools for code formatting and linting:

```bash
# Format code
black src/phaselab tests
isort src/phaselab tests

# Lint code
ruff check src/phaselab tests

# Type checking
mypy src/phaselab
```

Configuration is in `pyproject.toml`. Key standards:
- Line length: 100 characters
- Python 3.9+ syntax
- Type hints encouraged

### Committing Changes

Write clear, descriptive commit messages:

```
Add coherence calculation for multi-guide panels

- Implement batch coherence computation
- Add test cases for edge conditions
- Update documentation
```

### Submitting a Pull Request

1. Push your branch to your fork:
```bash
git push origin feature/your-feature-name
```

2. Open a Pull Request on GitHub

3. Fill out the PR template:
   - Describe what the PR does
   - Reference any related issues
   - Note any breaking changes

4. Wait for CI checks to pass

5. Address any review feedback

## Code Guidelines

### Module Structure

```
src/phaselab/
    __init__.py          # Public API exports
    core/                # Core utilities and constants
    crispr/              # CRISPR guide design
    quantum/             # Quantum coherence computation
    circadian/           # Circadian rhythm analysis
    therapy/             # Therapeutic optimization
```

### Adding New Features

1. **Core functionality** goes in the appropriate submodule
2. **Public API** should be exported in `__init__.py`
3. **Tests** are required for all new code
4. **Documentation** should be updated

### Testing Standards

- Aim for >80% code coverage
- Test edge cases and error conditions
- Use fixtures for common test data
- Mock external dependencies (quantum hardware, APIs)

Example test structure:
```python
import pytest
from phaselab import go_no_go

class TestGoNoGo:
    def test_above_threshold_is_go(self):
        result = go_no_go(0.5)
        assert result.status == "GO"

    def test_below_threshold_is_nogo(self):
        result = go_no_go(0.1)
        assert result.status == "NO-GO"

    def test_at_threshold_is_go(self):
        from phaselab.core.constants import E_MINUS_2
        result = go_no_go(E_MINUS_2)
        assert result.status == "GO"
```

### Documentation Standards

- All public functions need docstrings
- Use NumPy-style docstrings
- Include examples where helpful

Example:
```python
def compute_coherence(phases: np.ndarray) -> float:
    """
    Compute coherence from phase array.

    Parameters
    ----------
    phases : np.ndarray
        Array of phases in radians, shape (N,)

    Returns
    -------
    float
        Coherence value R-bar in range [0, 1]

    Examples
    --------
    >>> phases = np.array([0.1, 0.2, 0.15])
    >>> compute_coherence(phases)
    0.9876
    """
```

## Types of Contributions

### Bug Reports

Open an issue with:
- Clear description of the bug
- Steps to reproduce
- Expected vs actual behavior
- Python version and OS
- PhaseLab version

### Feature Requests

Open an issue describing:
- The use case
- Proposed solution
- Alternatives considered

### Documentation

Improvements welcome:
- Fix typos
- Clarify explanations
- Add examples
- Improve tutorials

### Code Contributions

- Bug fixes
- New features
- Performance improvements
- Test coverage improvements

## Optional Dependencies

When adding features that require optional dependencies:

1. Add to the appropriate group in `pyproject.toml`:
```toml
[project.optional-dependencies]
quantum = ["qiskit>=1.0", ...]
```

2. Use conditional imports:
```python
try:
    from qiskit import QuantumCircuit
    HAS_QISKIT = True
except ImportError:
    HAS_QISKIT = False
```

3. Skip tests when dependencies unavailable:
```python
@pytest.mark.skipif(not HAS_QISKIT, reason="qiskit not installed")
def test_quantum_feature():
    ...
```

## Release Process

Releases are managed by maintainers:

1. Version bump in `pyproject.toml`
2. Update CHANGELOG.md
3. Create git tag
4. GitHub Actions handles PyPI publishing

## Getting Help

- Open an issue for questions
- Check existing issues and discussions
- Read the documentation at https://phaselab.readthedocs.io

## Code of Conduct

Be respectful and constructive. We're all here to advance science.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
