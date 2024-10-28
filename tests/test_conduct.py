from pathlib import Path

from philharmonic.conduct import (
    build_snakemake_command,
    download_snakefile,
    get_snakefile_remote_path,
)


# Test get_snakefile_remote_path
def test_get_snakefile_remote_path():
    path = get_snakefile_remote_path()
    assert "https://raw.githubusercontent.com/samsledje/philharmonic/" in path
    assert path.endswith("Snakefile")

    path_with_commit = get_snakefile_remote_path("abc123")
    assert "abc123" in path_with_commit


# Test download_snakefile (this test might need to be adjusted based on your actual implementation)
def test_download_snakefile(tmp_path, monkeypatch):
    def mock_download_file_safe(url, output_path):
        output_path.write_text("mock snakefile content")
        return True

    monkeypatch.setattr(
        "philharmonic.conduct.download_file_safe", mock_download_file_safe
    )

    download_loc = tmp_path / "Snakefile"
    result = download_snakefile(download_loc)

    assert result == download_loc
    assert download_loc.exists()
    assert download_loc.read_text() == "mock snakefile content"


# Test build_snakemake_command
def test_build_snakemake_command():
    snakefile = Path("/path/to/Snakefile")
    config = Path("/path/to/config.yml")
    cores = 4

    cmd = build_snakemake_command(snakefile, config, cores)

    assert cmd[0] == "snakemake"
    assert "-s" in cmd and str(snakefile) in cmd
    assert "--configfile" in cmd and str(config) in cmd
    assert "--cores" in cmd and "4" in cmd
