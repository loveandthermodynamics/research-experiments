# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Verify delivered file bytes; this does not validate scientific conclusions."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    manifest = json.loads((ROOT / "release-integrity.json").read_text())
    for entry in manifest["files"]:
        relative = Path(entry["path"])
        path = (ROOT / relative).resolve()
        if relative.is_absolute() or not path.is_relative_to(ROOT):
            raise AssertionError("Integrity manifest contains a path outside the archive")
        if not path.is_file() or path.stat().st_size != entry["bytes"]:
            raise AssertionError(f"Missing file or changed size: {relative}")
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise AssertionError(f"Changed file bytes: {relative}")
    print(json.dumps({
        "files_checked": len(manifest["files"]),
        "all_matching": True,
        "scope": "File integrity only; not reproduction or scientific validation.",
    }, indent=2))


if __name__ == "__main__":
    main()
