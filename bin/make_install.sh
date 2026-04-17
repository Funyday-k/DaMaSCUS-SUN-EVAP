#!/bin/bash -l
cd ../build
cmake --build . --config Release
cmake --install .
cd ../bin