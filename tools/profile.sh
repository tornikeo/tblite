#!/bin/bash
set -e
meson compile -C build --verbose
timestamp=$(date +%Y%m%d_%H%M%S)
cd build && ncu --set detailed -f -o ../notes/profiles/$timestamp --open-in-ui test/unit/tester hamiltonian