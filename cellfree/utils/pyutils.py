from pathlib import Path


def examine_path(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return 0
