from pathlib import Path


def directory_path(path: str) -> Path:
    p = Path(path)
    if p.is_dir():
        return p
    raise NotADirectoryError(path)
