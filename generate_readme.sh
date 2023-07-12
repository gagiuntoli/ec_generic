#!/bin/bash

echo '# ec_generic'$'\n' > README.md
sed -n '/^\/\*!/,/\*\//p' src/elliptic_curve.rs >> README.md
sed -i '3d;$d' README.md
