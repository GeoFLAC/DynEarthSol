#!/bin/sh
set -eu

des_test_exe=$1
"$des_test_exe"
# A run with option $1 must stop with status 54, naming element 0 and the cause $2.
expect_rejection() {
    des_rejection_log="${des_test_exe}.rejected.log"
    if "$des_test_exe" "$1" > "$des_rejection_log" 2>&1; then
        cat "$des_rejection_log"
        echo "FAIL: a rejected Mohr-Coulomb return continued successfully ($1)" >&2
        exit 1
    else
        des_rejection_status=$?
    fi
    if [ "$des_rejection_status" -ne 54 ] ||
       ! grep -q 'no admissible Mohr-Coulomb return' "$des_rejection_log" ||
       ! grep -q "First element 0: $2" "$des_rejection_log"; then
        cat "$des_rejection_log"
        echo "FAIL: expected the cause '$2', element 0 and exit status 54 ($1)" >&2
        exit 1
    fi
    echo "PASS: a rejected Mohr-Coulomb return reports its cause and exit status 54 ($1)"
}
# In 3D the pure-2D option runs the same 3D case.
for des_option in --reject-mc-return --reject-mc-return-pure-2d; do
    expect_rejection "$des_option" 'non-finite trial stress'
done
# A law that jumps up where the return lands: no increment records itself.
for des_option in --reject-mc-jump --reject-mc-cap-jump; do
    expect_rejection "$des_option" 'no plastic increment returns the stress'
done
