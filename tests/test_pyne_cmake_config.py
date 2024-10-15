import subprocess
import pytest

def test_cmake_find_pyne(tmp_path):
    """
    Test if CMake can find the PyNE package correctly.
    """
    # Define the path to the temporary build directory.
    build_dir = tmp_path / "build"
    build_dir.mkdir()

    # Create a minimal CMakeLists.txt file in the temporary directory.
    cmake_content = """
    cmake_minimum_required(VERSION 3.15)
    project(TestPyNE)

    # Find PyNE package
    find_package(PyNE REQUIRED)
    """
    cmake_file = tmp_path / "CMakeLists.txt"
    cmake_file.write_text(cmake_content)

    # Run CMake in the temporary build directory.
    try:
        result = subprocess.run(
            ["cmake", str(tmp_path)],
            cwd=build_dir,
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        pytest.fail(f"CMake configuration failed with error: {e.stderr}")

    # Debug output
    print("CMake output:\n", result.stdout)
    pytest.warns(UserWarning, result.stderr)
