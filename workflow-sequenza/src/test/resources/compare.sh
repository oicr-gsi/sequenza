#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail
set -o noclobber

diff -s <(sort "$1") <(sort "$2")
