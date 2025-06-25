import unittest
import os
import tempfile
import shutil
import subprocess
from unittest.mock import patch

from .. import amalgamate


class TestVersionHandling(unittest.TestCase):
    """
    Tests the logic for extracting the version string from
    Git or archival files.
    """

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        amalgamate.BASE_DIR = self.test_dir

    def tearDown(self):
        # Clean up the directory
        shutil.rmtree(self.test_dir)

    @patch("subprocess.check_output")
    def test_get_version_from_git(self, mock_check_output):
        # Simulate successful git describe
        mock_check_output.side_effect = [
            b"true\n",  # is-inside-work-tree
            b"v1.2.3-10-gabc123\n",  # describe --tags
        ]
        version = amalgamate.get_version()
        self.assertEqual(version, "v1.2.3-10-gabc123")
        self.assertEqual(mock_check_output.call_count, 2)

    @patch("subprocess.check_output")
    def test_get_version_git_fails(self, mock_check_output):
        # Simulate git describe failing
        mock_check_output.side_effect = [
            b"true\n",  # is-inside-work-tree
            subprocess.CalledProcessError(1, "git", output=b"fatal: No names found"),
        ]
        with self.assertRaisesRegex(RuntimeError, "Git describe failed"):
            amalgamate.get_version()

    @patch("subprocess.check_output", side_effect=FileNotFoundError)
    def test_get_version_from_archival(self, mock_check_output):
        # Simulate not being in a git repo, but having an archival file
        archival_path = os.path.join(self.test_dir, ".git_archival.txt")
        with open(archival_path, "w") as f:
            f.write("ref-names: tags/v1.2.3\n")
            f.write("describe-name: v1.2.3-archived\n")

        version = amalgamate.get_version()
        self.assertEqual(version, "v1.2.3-archived")

    @patch("subprocess.check_output", side_effect=FileNotFoundError)
    def test_get_version_no_git_no_archival(self, mock_check_output):
        # Simulate worst case: no git, no archival file
        with self.assertRaisesRegex(
            RuntimeError, "Not in a Git repo and .git_archival.txt is missing."
        ):
            amalgamate.get_version()

    def test_create_version_header(self):
        test_version = "1.2.3-test"
        test_filename = "custom_version_output.h"

        # Calculate the full path where we expect the file to be created.
        expected_path = os.path.join(self.test_dir, test_filename)

        # Temporarily patch `amalgamate.BASE_DIR` to point to test directory.
        with patch.object(amalgamate, "BASE_DIR", self.test_dir):
            amalgamate.create_version_header(test_version, test_filename)

        # Assert that the file was actually created at the expected path.
        self.assertTrue(
            os.path.exists(expected_path),
            f"File was not created at expected path: {expected_path}",
        )

        # Read the created file and verify its content.
        with open(expected_path, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertIn("#ifndef PYNE_VERSION_HEADER", content)
        self.assertIn(f'return "{test_version} (amalgamated)";', content)


class TestAmalgamatedFile(unittest.TestCase):
    """
    Tests the core class responsible for assembling the amalgamated
    files, including file concatenation and include filtering.
    """

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        # Create some dummy files
        self.header_file = os.path.join(self.test_dir, "test.h")
        with open(self.header_file, "w") as f:
            f.write("void test_func();\n")

        self.source_file = os.path.join(self.test_dir, "test.cpp")
        with open(self.source_file, "w") as f:
            f.write('#include "test.h"\n')

        self.text_file = os.path.join(self.test_dir, "license.txt")
        with open(self.text_file, "w") as f:
            f.write("This is a license.\n")

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_append_file_and_write(self):
        output_path = os.path.join(self.test_dir, "output.cpp")
        am_file = amalgamate.AmalgamatedFile(output_path)
        am_file.append_file(self.header_file)
        am_file.append_file(self.source_file)
        am_file.write()

        self.assertTrue(os.path.exists(output_path))
        with open(output_path, "r") as f:
            content = f.read()

        self.assertIn(f"// Begin: {self.header_file}", content)
        self.assertIn("void test_func();", content)
        self.assertIn(f"// End: {self.header_file}", content)
        self.assertIn(f"// Begin: {self.source_file}", content)
        self.assertIn('#include "test.h"', content)
        self.assertIn(f"// End: {self.source_file}", content)

    def test_non_code_file_is_commented(self):
        output_path = os.path.join(self.test_dir, "output.h")
        am_file = amalgamate.AmalgamatedFile(output_path)
        am_file.append_file(self.text_file)
        am_file.write()

        with open(output_path, "r") as f:
            content = f.read()

        self.assertIn("// This is a license.", content)
        # Check that the line isn't present without the comment
        self.assertNotIn("\nThis is a license.", content)

    def test_prepend_file_listing(self):
        output_path = os.path.join(self.test_dir, "output.cpp")
        am_file = amalgamate.AmalgamatedFile(output_path)
        am_file.append_file(self.source_file)
        am_file.append_file(self.header_file)
        am_file.write()

        with open(output_path, "r") as f:
            content = f.read()

        self.assertTrue(
            content.startswith("// Amalgamated from the following files:\n")
        )
        self.assertIn(f"//   {self.source_file}\n", content)
        self.assertIn(f"//   {self.header_file}\n", content)

    def test_write_creates_directory(self):
        output_dir = os.path.join(self.test_dir, "new_dir")
        output_path = os.path.join(output_dir, "output.cpp")

        self.assertFalse(os.path.exists(output_dir))

        am_file = amalgamate.AmalgamatedFile(output_path)
        am_file.append_line("int main() { return 0; }")
        am_file.write()

        self.assertTrue(os.path.exists(output_dir))
        self.assertTrue(os.path.exists(output_path))


class TestMainExecution(unittest.TestCase):
    """
    Tests the main script logic, argument parsing, and
    overall file generation process using mock files.
    """

    def setUp(self):
        self.original_cwd = os.getcwd()
        self.test_dir = tempfile.mkdtemp()
        # Change to the test directory to handle relative paths like "version.h" cleanly
        os.chdir(self.test_dir)

        self.src_dir = os.path.join(self.test_dir, "src")
        os.makedirs(self.src_dir)

        # Create dummy source files
        self.files_map = {
            "license.txt": "PyNE License",
            "src/utils.h": "namespace pyne { int util(); }",
            "src/utils.cpp": '#include "utils.h"\nint pyne::util() { return 42; }',
            "src/nucname.h": "namespace pyne { void nucname(); }",
            "src/nucname.cpp": '#include "nucname.h"\nvoid pyne::nucname() {}',
        }
        self.file_paths = []
        for path, content in self.files_map.items():
            full_path = os.path.join(self.test_dir, path)
            os.makedirs(os.path.dirname(full_path), exist_ok=True)
            with open(full_path, "w") as f:
                f.write(content)
            self.file_paths.append(full_path)

        # The main() function creates version.h, so we must add it to the list of files to process.
        self.file_paths.append("version.h")

    def tearDown(self):
        # Change back to the original directory before removing the temp one
        os.chdir(self.original_cwd)
        shutil.rmtree(self.test_dir)

    @patch.object(amalgamate, "get_version", return_value="9.9.9-test")
    @patch("sys.argv", new_callable=list)
    def test_main_default_args(self, mock_argv, mock_get_version):
        # Setup mocks and arguments
        amalgamate.BASE_DIR = self.test_dir

        # Note: output dir is now '.' because we chdir'd into self.test_dir
        mock_argv.extend(["amalgamate.py", "-o", ".", "-f", *self.file_paths])

        # Run the main function
        amalgamate.main()

        # Assertions
        header_path = os.path.join(self.test_dir, "pyne.h")
        source_path = os.path.join(self.test_dir, "pyne.cpp")

        # Check files were created
        self.assertTrue(os.path.exists(header_path))
        self.assertTrue(os.path.exists(source_path))

        # Check temporary version.h was deleted
        self.assertFalse(os.path.exists("version.h"))

        # Check header content
        with open(header_path, "r") as f:
            header_content = f.read()
        self.assertIn("#ifndef PYNE_AMALGAMATED_HEADER", header_content)
        self.assertIn("#define PYNE_IS_AMALGAMATED", header_content)
        self.assertIn(self.files_map["src/utils.h"], header_content)
        self.assertIn(self.files_map["src/nucname.h"], header_content)
        self.assertNotIn(
            self.files_map["src/utils.cpp"], header_content
        )  # Cpp not in header
        self.assertIn("// " + self.files_map["license.txt"], header_content)
        self.assertIn('return "9.9.9-test (amalgamated)";', header_content)

        # Check source content
        with open(source_path, "r") as f:
            source_content = f.read()
        self.assertIn('#include "pyne.h"', source_content)
        self.assertIn(self.files_map["src/utils.cpp"], source_content)
        self.assertIn(self.files_map["src/nucname.cpp"], source_content)
        self.assertNotIn(
            self.files_map["src/utils.h"], source_content
        )  # Header not in source
        # This assertion will now pass
        self.assertIn("// " + self.files_map["license.txt"], source_content)


class TestAmalgamationResult(unittest.TestCase):
    """
    An end-to-end integration test that checks if the amalgamated
    output files can be compiled and used in a real C++ program.
    """

    def setUp(self):
        # Create a top-level temporary directory for the test
        self.test_dir = tempfile.mkdtemp()
        self.output_dir = os.path.join(self.test_dir, "amalgamated_lib")
        os.makedirs(self.output_dir)

        # Patch BASE_DIR for the amalgamate script
        self.base_dir_patcher = patch.object(amalgamate, "BASE_DIR", self.test_dir)
        self.base_dir_patcher.start()

    def tearDown(self):
        self.base_dir_patcher.stop()
        shutil.rmtree(self.test_dir)

    @unittest.skipIf(
        shutil.which("cmake") is None,
        "CMake executable not found, skipping build test.",
    )
    def test_build_and_run_amalgamated_library(self):
        """
        Tests the full cycle: create sources -> amalgamate -> compile -> run.
        """
        # Create dummy source files for the library
        src_dir = os.path.join(self.test_dir, "src")
        os.makedirs(src_dir)

        lib_header_path = os.path.join(src_dir, "functions.h")
        lib_source_path = os.path.join(src_dir, "functions.cpp")
        version_header_path = os.path.join(
            self.test_dir, "version.h"
        )  # Will be created by main()

        with open(lib_header_path, "w") as f:
            f.write(
                "#include <string>\n" "namespace pyne { std::string get_greeting(); }"
            )

        with open(lib_source_path, "w") as f:
            f.write(
                '#include "functions.h"\n'
                'namespace pyne { std::string get_greeting() { return "Hello from amalgamated PyNE!"; } }'
            )

        # Run the amalgamation script via its main() function
        source_files_to_process = [
            lib_header_path,
            lib_source_path,
            version_header_path,
        ]

        with patch(
            "sys.argv",
            ["amalgamate.py", "-o", self.output_dir, "-f", *source_files_to_process],
        ):
            with patch.object(amalgamate, "get_version", return_value="2.0-cmake-test"):
                amalgamate.main()

        # Create CMake and a C++ test program to use the library
        main_test_cpp = os.path.join(self.output_dir, "main_test.cpp")
        cmakelists_txt = os.path.join(self.output_dir, "CMakeLists.txt")

        # This C++ program will fail if the library functions don't work
        with open(main_test_cpp, "w") as f:
            f.write(
                """
#include <iostream>
#include <string>
#include "pyne.h" // The amalgamated header

int main() {
    std::string greeting = pyne::get_greeting();
    std::string version = pyne::pyne_version();

    if (greeting != "Hello from amalgamated PyNE!") {
        std::cerr << "FAIL: get_greeting() returned incorrect value." << std::endl;
        return 1;
    }

    if (version.find("2.0-cmake-test") == std::string::npos) {
        std::cerr << "FAIL: pyne_version() returned incorrect value." << std::endl;
        return 1;
    }

    std::cout << "SUCCESS: Functional test passed." << std::endl;
    std::cout << "Greeting: " << greeting << std::endl;
    std::cout << "Version: " << version << std::endl;
    return 0;
}
"""
            )
        # This CMake file builds the test program
        with open(cmakelists_txt, "w") as f:
            f.write(
                """
cmake_minimum_required(VERSION 3.10)
project(AmalgamationFunctionTest CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add our test program, which is composed of the test runner
# and the amalgamated source file.
add_executable(run_functional_test main_test.cpp pyne.cpp)
"""
            )

        # Run CMake, build, and execute the test program
        build_dir = os.path.join(self.output_dir, "build")
        os.makedirs(build_dir)

        try:
            # Configure with CMake
            subprocess.run(
                ["cmake", ".."],
                cwd=build_dir,
                check=True,
                capture_output=True,
                text=True,
            )
            # Build with CMake
            subprocess.run(
                ["cmake", "--build", "."],
                cwd=build_dir,
                check=True,
                capture_output=True,
                text=True,
            )
            # Run the compiled executable
            executable_name = "run_functional_test"
            if os.name == "nt":
                executable_name += ".exe"

            result = subprocess.run(
                [os.path.join(build_dir, executable_name)],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            # Provide detailed error output if any command fails
            self.fail(
                f"CMake build or run failed.\n"
                f"STDOUT:\n{e.stdout}\n"
                f"STDERR:\n{e.stderr}"
            )

        # Assert that the output from the C++ program is correct
        print("Executable Output:\n---", result.stdout, "---")
        self.assertIn("SUCCESS: Functional test passed.", result.stdout)
        self.assertIn("Hello from amalgamated PyNE!", result.stdout)
        self.assertIn("2.0-cmake-test (amalgamated)", result.stdout)


if __name__ == "__main__":
    unittest.main(argv=["first-arg-is-ignored"], exit=False)
