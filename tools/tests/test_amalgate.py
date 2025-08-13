import pytest
import sys
import subprocess
from pathlib import Path
from unittest.mock import patch

# Import the script to be tested.
from .. import amalgamate

@pytest.fixture
def test_env(tmp_path):
    """
    Creates a temporary file system structure for tests to run in isolation.
    This fixture provides a realistic environment for the amalgamation script.
    """
    # Create base directory for the test project
    base_dir = tmp_path / "pyne_project"
    base_dir.mkdir()

    # Create dummy source and header files
    src_dir = base_dir / "src"
    src_dir.mkdir()
    (src_dir / "utils.h").write_text('// Content of utils.h\n#include <string>\nstd::string get_util_name();\n')
    (src_dir / "utils.cpp").write_text('// Content of utils.cpp\n#include "utils.h"\nstd::string get_util_name() { return "util"; }\n')
    (src_dir / "data.h").write_text('// Content of data.h\nconst int THE_ANSWER = 42;\n')
    (src_dir / "data.cpp").write_text('// Content of data.cpp\n#include "data.h"\n')
    (src_dir / "unsupported.txt").write_text('This file should be skipped.')

    # Create a dummy license file
    (base_dir / "license.txt").write_text("The PyNE License Block.")

    # Create a git archival file for versioning fallback
    (base_dir / ".git_archival.txt").write_text("describe-name: v0.7.7-archival")

    # The amalgamate script itself is considered to be at the base level
    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(amalgamate, "BASE_DIR", base_dir)

    yield {
        "base": base_dir,
        "src": src_dir,
        "output": tmp_path / "output",
    }
    
    monkeypatch.undo()

# Test get_version function
def test_get_version_from_git(monkeypatch):
    """Tests that get_version() correctly calls git and returns its output."""
    # Mock subprocess.check_output to simulate git commands
    def mock_git(cmd, **kwargs):
        if "rev-parse" in cmd:
            return b'true\n'
        if "describe" in cmd:
            return b'v1.2.3-git\n'
        return b''

    monkeypatch.setattr(subprocess, "check_output", mock_git)
    version = amalgamate.get_version()
    assert version == "v1.2.3-git"

def test_get_version_from_archival(test_env, monkeypatch):
    """Tests that get_version() falls back to .git_archival.txt when git is not available."""
    # Mock subprocess.check_output to fail, simulating a non-git environment
    def mock_fail(*args, **kwargs):
        raise FileNotFoundError("git not found")

    monkeypatch.setattr(subprocess, "check_output", mock_fail)
    version = amalgamate.get_version()
    assert version == "v0.7.7-archival"

def test_get_version_git_describe_fails(monkeypatch):
    """Tests that get_version() raises RuntimeError if 'git describe' fails."""
    def mock_git_fail(cmd, **kwargs):
        if "rev-parse" in cmd:
            return b'true\n'
        if "describe" in cmd:
            raise subprocess.CalledProcessError(1, cmd, output=b"fatal: No tags can describe 'HEAD'.")
        return b''

    monkeypatch.setattr(subprocess, "check_output", mock_git_fail)
    with pytest.raises(RuntimeError) as excinfo:
        amalgamate.get_version()
    assert "Git describe failed" in str(excinfo.value)

def test_get_version_fails_without_git_or_archival(test_env, monkeypatch):
    """Tests that get_version() fails if there's no version source."""
    (test_env["base"] / ".git_archival.txt").unlink() # Remove the fallback file
    
    def mock_fail(*args, **kwargs):
        raise FileNotFoundError("git not found")

    monkeypatch.setattr(subprocess, "check_output", mock_fail)
    with pytest.raises(RuntimeError) as excinfo:
        amalgamate.get_version()
    assert "Not in a Git repo and .git_archival.txt is missing" in str(excinfo.value)

# Test main script execution
def test_full_run_with_defaults(test_env, monkeypatch):
    """Tests a standard execution of the script, checking file content and structure."""
    # Setup command-line arguments using monkeypatch
    input_files = [
        test_env["src"] / "utils.h",
        test_env["src"] / "utils.cpp",
        test_env["src"] / "data.h",
        test_env["src"] / "data.cpp",
    ]
    args = [
        "amalgamate.py",
        "-o", str(test_env["output"]),
        "-f", *[str(f) for f in input_files],
    ]
    monkeypatch.setattr(sys, "argv", args)
    monkeypatch.setattr(amalgamate, "get_version", lambda: "v-test")

    # Run the main function
    amalgamate.main()

    # Verification
    header_file = test_env["output"] / "pyne.h"
    source_file = test_env["output"] / "pyne.cpp"

    assert header_file.exists()
    assert source_file.exists()

    header_content = header_file.read_text()
    source_content = source_file.read_text()

    # Check for license and version
    assert "The PyNE License Block." in header_content
    assert 'return "v-test (amalgamated)";' in header_content
    assert "The PyNE License Block." in source_content

    # Check for include guards and definition
    assert "#ifndef PYNE_AMALGAMATED_HEADER" in header_content
    assert "#define PYNE_IS_AMALGAMATED" in header_content
    assert "#endif  // PYNE_AMALGAMATED_HEADER" in header_content

    # Check for correct file inclusion
    assert "// Content of utils.h" in header_content
    assert "// Content of data.h" in header_content
    assert "// Content of utils.cpp" not in header_content

    assert "// Content of utils.cpp" in source_content
    assert "// Content of data.cpp" in source_content
    assert "// Content of data.h" not in source_content
    
    # Check that source file includes the header
    assert '#include "pyne.h"' in source_content

    # Check for file listing comment
    assert "Amalgamated from the following files:" in header_content
    assert "src/utils.h" in header_content
    assert "Amalgamated from the following files:" in source_content
    assert "src/utils.cpp" in source_content

def test_custom_output_filenames(test_env, monkeypatch):
    """Tests the -s (source) and -i (header) arguments for custom file naming."""
    args = [
        "amalgamate.py",
        "-o", str(test_env["output"]),
        "-s", "custom_pyne.cc",
        "-i", "custom_pyne.hh",
        "-f", str(test_env["src"] / "utils.h"),
    ]
    monkeypatch.setattr(sys, "argv", args)
    monkeypatch.setattr(amalgamate, "get_version", lambda: "v-custom")

    amalgamate.main()

    custom_header = test_env["output"] / "custom_pyne.hh"
    custom_source = test_env["output"] / "custom_pyne.cc"

    assert custom_header.exists()
    assert custom_source.exists()
    
    # Check that the source includes the custom-named header
    source_content = custom_source.read_text()
    assert '#include "custom_pyne.hh"' in source_content

def test_relative_header_include_path(test_env, monkeypatch):
    """Tests that the #include path is correctly relative when outputs are in different directories."""
    # Define separate, nested output paths for header and source
    header_path = test_env["output"] / "include" / "pyne.h"
    source_path = test_env["output"] / "source" / "pyne.cpp"
    
    args = [
        "amalgamate.py",
        "-o", ".", # Output relative to current dir
        "-i", str(header_path),
        "-s", str(source_path),
        "-f", str(test_env["src"] / "data.cpp"),
    ]
    monkeypatch.setattr(sys, "argv", args)
    monkeypatch.setattr(amalgamate, "get_version", lambda: "v-relative")

    amalgamate.main()
    
    assert header_path.exists()
    assert source_path.exists()
    
    source_content = source_path.read_text()
    # Path from output/source/pyne.cpp to output/include/pyne.h is ../include/pyne.h
    assert '#include "../include/pyne.h"' in source_content

def test_missing_license_warning(test_env, monkeypatch, capsys):
    """Tests that a warning is printed and handled if license.txt is missing."""
    (test_env["base"] / "license.txt").unlink() # Delete the license file

    monkeypatch.setattr(sys, "argv", ["amalgamate.py", "-o", str(test_env["output"])])
    monkeypatch.setattr(amalgamate, "get_version", lambda: "v-nolicense")

    amalgamate.main()

    # Check for the warning message in stdout
    captured = capsys.readouterr()
    assert "[!] Warning: license.txt not found." in captured.out

    # Check that the output files contain the placeholder comment
    header_content = (test_env["output"] / "pyne.h").read_text()
    assert "// license.txt not found." in header_content

def test_skip_unknown_extension_warning(test_env, monkeypatch, capsys):
    """Tests that a warning is printed for files with unsupported extensions."""
    args = [
        "amalgamate.py",
        "-o", str(test_env["output"]),
        "-f", str(test_env["src"] / "unsupported.txt"),
    ]
    monkeypatch.setattr(sys, "argv", args)
    monkeypatch.setattr(amalgamate, "get_version", lambda: "v-skip")

    amalgamate.main()

    # Check for the warning message in stdout
    captured = capsys.readouterr()
    assert "[!] Warning: Skipping file with unknown extension: " in captured.out
    assert "unsupported.txt" in captured.out