# Note: when executed in the build dir, then CMAKE_CURRENT_SOURCE_DIR is the
# build dir.
file(COPY scripts DESTINATION "${CMAKE_ARGV3}")
file(COPY setup.py setup_sub.py pyne DESTINATION "${CMAKE_ARGV3}"
     FILES_MATCHING PATTERN "*.py"
                    PATTERN "*.pyw"
                    PATTERN "*.csv"
                    PATTERN "*.txt"
                    PATTERN "*.inp"
                    PATTERN "*.html"
                    PATTERN "*.pxi"
                    PATTERN "*.pyi"
                    PATTERN "*.pxd"
                    PATTERN "*.mas16"
                    PATTERN "*.dat"
                    )
