// Copyright 2006, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <iostream>

#include "gtest/gtest.h"

char ** gtest_argv;
int gtest_argc;

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running main() from gtest_main.cc\n";

  // jcz -- Provide the capability of passing arguments
  // ref -- groups.google.com/forum
  const char * data_path = NULL;
  gtest_argv = (char **)malloc( sizeof(char *) * argc);
  gtest_argv[0] = argv[0]; // copying the first argument
  gtest_argc = 1;
  for (int i=1; i<argc; i++)
  {
      if (strcmp(argv[i],"data_path") != 0)
      {
	 gtest_argv[gtest_argc++] = argv[i];
      }
      else
      {
         data_path = argv[i];
      }
  }
  ///////////////// end of additional code

  // testing::InitGoogleTest(&argc, argv);
  testing::InitGoogleTest(&gtest_argc, gtest_argv);
  free(gtest_argv);
  return RUN_ALL_TESTS();
}
