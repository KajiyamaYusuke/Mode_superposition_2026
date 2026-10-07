"""Canonical result layout with read compatibility for legacy flat outputs."""

from pathlib import Path

EXTENSION_DIR = {
    ".csv": "csv", ".dat": "dat", ".png": "png", ".txt": "txt",
    ".wav": "wav", ".json": "json", ".pdf": "pdf",
}


def result_path(run_dir: Path, name: str | Path) -> Path:
    """Prefer the organized path, but continue to read legacy flat outputs."""
    name = Path(name)
    subdir = EXTENSION_DIR.get(name.suffix.lower())
    canonical = run_dir / subdir / name if subdir else run_dir / name
    legacy = run_dir / name
    return canonical if canonical.is_file() or not legacy.is_file() else legacy


def output_path(run_dir: Path, name: str | Path) -> Path:
    """Return the organized output path and create its parent directory."""
    name = Path(name)
    subdir = EXTENSION_DIR.get(name.suffix.lower())
    path = run_dir / subdir / name if subdir else run_dir / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
