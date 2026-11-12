#!/bin/sh
set -eu

des_test_exe=$1
"$des_test_exe"
des_rejection_log="${des_test_exe}.rejected.log"
if "$des_test_exe" --reject-mc-return > "$des_rejection_log" 2>&1; then
    cat "$des_rejection_log"
    echo 'FAIL: rejected Mohr-Coulomb return continued successfully' >&2
    exit 1
else
    des_rejection_status=$?
fi
if [ "$des_rejection_status" -ne 54 ] ||
   ! grep -q 'no admissible Mohr-Coulomb return' "$des_rejection_log" ||
   ! grep -q 'First element 0: non-finite trial stress' "$des_rejection_log"; then
    cat "$des_rejection_log"
    echo 'FAIL: expected the non-finite cause, element 0 and exit status 54' >&2
    exit 1
fi
echo 'PASS: rejected Mohr-Coulomb return reports its cause and exit status 54'
