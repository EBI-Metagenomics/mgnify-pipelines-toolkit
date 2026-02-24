import io
import json
from pathlib import Path
from urllib.error import HTTPError, URLError

import pytest

from mgnify_pipelines_toolkit.analysis.shared.download_zenodo import (
    download_file,
    get_zenodo_metadata,
    main,
)


class _FakeHTTPResponse:
    """Context-manager response object compatible with urllib.request.urlopen()."""

    def __init__(self, payload: bytes):
        self._payload = payload

    def read(self) -> bytes:
        return self._payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def test_download_file_success(monkeypatch, tmp_path: Path):
    # Arrange: mock urlopen to return bytes
    def fake_urlopen(_req, timeout=None):
        return _FakeHTTPResponse(b"hello-world")

    monkeypatch.setattr("mgnify_pipelines_toolkit.analysis.shared.download_zenodo.urlopen", fake_urlopen)

    out = tmp_path / "file.bin"

    # Act
    ok = download_file(
        url="https://example.org/file.bin",
        output_file=str(out),
        user_agent="pytest-agent",
        max_retries=3,
        retry_delay=0,
    )

    # Assert
    assert ok is True
    assert out.read_bytes() == b"hello-world"


def test_download_file_retries_then_success(monkeypatch, tmp_path: Path):
    calls = {"n": 0}

    def fake_urlopen(_req, timeout=None):
        calls["n"] += 1
        if calls["n"] < 3:
            raise URLError("temporary failure")
        return _FakeHTTPResponse(b"ok-after-retry")

    monkeypatch.setattr("mgnify_pipelines_toolkit.analysis.shared.download_zenodo.urlopen", fake_urlopen)

    # Avoid real sleep
    monkeypatch.setattr(
        "mgnify_pipelines_toolkit.analysis.shared.download_zenodo.time.sleep",
        lambda _s: None,
    )

    out = tmp_path / "file.bin"

    ok = download_file(
        url="https://example.org/file.bin",
        output_file=str(out),
        user_agent="pytest-agent",
        max_retries=5,
        retry_delay=0,
    )

    assert ok is True
    assert calls["n"] == 3
    assert out.read_bytes() == b"ok-after-retry"


def test_download_file_exhausts_retries(monkeypatch, tmp_path: Path, capsys):
    def fake_urlopen(_req, timeout=None):
        raise URLError("permanent failure")

    monkeypatch.setattr("mgnify_pipelines_toolkit.analysis.shared.download_zenodo.urlopen", fake_urlopen)
    monkeypatch.setattr(
        "mgnify_pipelines_toolkit.analysis.shared.download_zenodo.time.sleep",
        lambda _s: None,
    )

    out = tmp_path / "file.bin"

    ok = download_file(
        url="https://example.org/file.bin",
        output_file=str(out),
        user_agent="pytest-agent",
        max_retries=3,
        retry_delay=0,
    )

    assert ok is False
    assert not out.exists()

    err = capsys.readouterr().err
    assert "Error downloading file after 3 attempts" in err


def test_get_zenodo_metadata_success(monkeypatch):
    payload = {
        "metadata": {"version": "1.2.3"},
        "files": [{"key": "db.tar.gz", "size": 123}],
    }

    def fake_urlopen(_req, timeout=None):
        return _FakeHTTPResponse(json.dumps(payload).encode("utf-8"))

    monkeypatch.setattr("mgnify_pipelines_toolkit.analysis.shared.download_zenodo.urlopen", fake_urlopen)

    data = get_zenodo_metadata(zenodo_id=12345, user_agent="pytest-agent")
    assert data["metadata"]["version"] == "1.2.3"
    assert data["files"][0]["key"] == "db.tar.gz"


def test_get_zenodo_metadata_http_error_exits(monkeypatch, capsys):
    # HTTPError(url, code, msg, hdrs, fp)
    def fake_urlopen(_req, timeout=None):
        raise HTTPError(
            url="https://zenodo.org/api/records/123",
            code=404,
            msg="Not Found",
            hdrs=None,
            fp=io.BytesIO(b""),
        )

    monkeypatch.setattr("mgnify_pipelines_toolkit.analysis.shared.download_zenodo.urlopen", fake_urlopen)

    with pytest.raises(SystemExit) as e:
        get_zenodo_metadata(zenodo_id=123, user_agent="pytest-agent")
    assert e.value.code == 1

    err = capsys.readouterr().err
    assert "HTTP Error fetching Zenodo metadata" in err


def test_main_happy_path(monkeypatch, tmp_path: Path, capsys):
    # Mock metadata + download_file, avoid any network and write to tmp_path
    payload = {
        "metadata": {"version": "9.9.9"},
        "files": [{"key": "db.tar.gz", "size": 10}],
    }

    def fake_get_meta(zenodo_id, user_agent):
        assert zenodo_id == 999
        return payload

    def fake_download(url, output_file, user_agent, max_retries=3, retry_delay=5):
        # emulate download by writing a small file
        Path(output_file).write_bytes(b"x" * 10)
        return True

    monkeypatch.setattr(
        "mgnify_pipelines_toolkit.analysis.shared.download_zenodo.get_zenodo_metadata",
        fake_get_meta,
    )
    monkeypatch.setattr(
        "mgnify_pipelines_toolkit.analysis.shared.download_zenodo.download_file",
        fake_download,
    )

    out = tmp_path / "out.tar.gz"
    monkeypatch.setattr(
        "mgnify_pipelines_toolkit.analysis.shared.download_zenodo.sys.argv",
        ["download_zenodo.py", "999", "--output-file", str(out), "--file-index", "0"],
    )

    # main() ends normally (no return), so just call it
    main()

    assert out.exists()
    assert out.read_bytes() == b"x" * 10

    stdout = capsys.readouterr().out
    assert "Fetching Zenodo record metadata for ID: 999" in stdout
    assert "Record version: 9.9.9" in stdout
    assert "Download complete" in stdout


def test_main_file_index_out_of_range_exits(monkeypatch, capsys):
    payload = {
        "metadata": {"version": "1.0.0"},
        "files": [{"key": "a.bin", "size": 1}],
    }

    monkeypatch.setattr(
        "mgnify_pipelines_toolkit.analysis.shared.download_zenodo.get_zenodo_metadata",
        lambda zenodo_id, user_agent: payload,
    )
    monkeypatch.setattr(
        "mgnify_pipelines_toolkit.analysis.shared.download_zenodo.sys.argv",
        ["download_zenodo.py", "123", "--file-index", "9"],
    )

    with pytest.raises(SystemExit) as e:
        main()
    assert e.value.code == 1

    err = capsys.readouterr().err
    assert "Error: File index 9 not found" in err
