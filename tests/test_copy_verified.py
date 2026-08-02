from pathlib import Path

import pytest

from workflow.scripts.copy_verified import copy_verified


def test_copy_verified_preserves_bytes(tmp_path):
    source = tmp_path / "source.bin"
    destination = tmp_path / "nested" / "destination.bin"
    source.write_bytes(bytes(range(256)) * 100)

    attempt = copy_verified(source, destination)

    assert attempt == 1
    assert destination.read_bytes() == source.read_bytes()


def test_copy_verified_replaces_existing_destination(tmp_path):
    source = tmp_path / "source.bin"
    destination = tmp_path / "destination.bin"
    source.write_bytes(b"expected")
    destination.write_bytes(b"stale")

    copy_verified(source, destination)

    assert destination.read_bytes() == b"expected"


def test_copy_verified_requires_positive_attempts(tmp_path):
    source = tmp_path / "source.bin"
    source.write_bytes(b"data")

    with pytest.raises(ValueError, match="attempts must be positive"):
        copy_verified(source, tmp_path / "destination.bin", attempts=0)
