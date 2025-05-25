#!/bin/bash
set -e
meson compile -C build --verbose
timestamp=$(date +%Y%m%d_%H%M%S)
# --set full
cd build && ncu -f -o ../notes/profiles/$timestamp -c 1 --open-in-ui test/unit/tester hamiltonian