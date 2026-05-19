#!/bin/bash
# Set up DES performance benchmark directories.
# Run from the benchmarks/DESv2.0/performance/ directory.
#
# Creates:
#   01-regular/    regular mesh throughput sweep  (CPU + GPU, 8 sizes)
#   02-irregular/  irregular mesh throughput sweep (CPU + GPU, 8 sizes)
#   03-scaling/    CPU core-scaling test (regular mesh, 4 sizes)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DES_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

REG_CFG="$SCRIPT_DIR/test-regular-3d.cfg"
IRREG_CFG="$SCRIPT_DIR/test-irregular-3d.cfg"

CPU_EXE="$DES_ROOT/dynearthsol3d"
GPU_EXE="$DES_ROOT/dynearthsol3d.gpu"

THROUGHPUT_SIZES=(2e3 4e3 8e3 16e3 32e3 64e3 128e3 256e3)
SCALING_SIZES=(2e3 4e3 8e3 16e3)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

require_file() {
    if [[ ! -f "$1" ]]; then
        echo "ERROR: required file not found: $1" >&2
        exit 1
    fi
}

gen_cfg() {
    local template="$1"
    local outfile="$2"
    local ylength="$3"
    local modelname="$4"
    sed -e "s|^ylength = .*|ylength = $ylength|" \
        -e "s|^modelname = .*|modelname = $modelname|" \
        "$template" > "$outfile"
}

copy_exes() {
    local dest="$1"
    local do_gpu="${2:-yes}"
    if [[ -f "$CPU_EXE" ]]; then
        cp "$CPU_EXE" "$dest/"
    else
        echo "  WARNING: CPU executable not found at $CPU_EXE"
    fi
    if [[ "$do_gpu" == "yes" ]]; then
        if [[ -f "$GPU_EXE" ]]; then
            cp "$GPU_EXE" "$dest/"
        else
            echo "  WARNING: GPU executable not found at $GPU_EXE"
        fi
    fi
}

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------

require_file "$REG_CFG"
require_file "$IRREG_CFG"

# ---------------------------------------------------------------------------
# 01-regular   (regular mesh, CPU + GPU throughput sweep)
# ---------------------------------------------------------------------------
DIR="$SCRIPT_DIR/01-regular"
echo "Creating $DIR ..."
mkdir -p "$DIR"

idx=0
for yl in "${THROUGHPUT_SIZES[@]}"; do
    idx=$((idx+1)); printf -v n "%02d" $idx
    gen_cfg "$REG_CFG" "$DIR/test-regular-3d-${n}-${yl}.cfg" "$yl" "reg-$yl"
done

copy_exes "$DIR" yes

cat > "$DIR/run-cpu.sh" << 'EOF'
#!/bin/bash
# Regular mesh — CPU throughput sweep
# Adjust OMP_NUM_THREADS to the number of physical cores on the test machine.
set -euo pipefail
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-64}
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
for cfg in test-regular-3d-*.cfg; do
    echo "--- $cfg ---"
    ./dynearthsol3d "$cfg" 2>&1 | tee "${cfg%.cfg}-cpu.log"
done
echo "Done. Collect: grep 's/step' *-cpu.log"
EOF

cat > "$DIR/run-gpu.sh" << 'EOF'
#!/bin/bash
# Regular mesh — GPU throughput sweep
set -euo pipefail
for cfg in test-regular-3d-*.cfg; do
    echo "--- $cfg ---"
    ./dynearthsol3d.gpu "$cfg" 2>&1 | tee "${cfg%.cfg}-gpu.log"
done
echo "Done. Collect: grep 's/step' *-gpu.log"
EOF

chmod +x "$DIR/run-cpu.sh" "$DIR/run-gpu.sh"
echo "  $(ls "$DIR"/*.cfg | wc -l) cfg files, run-cpu.sh + run-gpu.sh"

# ---------------------------------------------------------------------------
# 02-irregular  (irregular mesh, CPU + GPU throughput sweep)
# ---------------------------------------------------------------------------
DIR="$SCRIPT_DIR/02-irregular"
echo "Creating $DIR ..."
mkdir -p "$DIR"

idx=0
for yl in "${THROUGHPUT_SIZES[@]}"; do
    idx=$((idx+1)); printf -v n "%02d" $idx
    gen_cfg "$IRREG_CFG" "$DIR/test-irregular-3d-${n}-${yl}.cfg" "$yl" "irreg-$yl"
done

copy_exes "$DIR" yes

cat > "$DIR/run-cpu.sh" << 'EOF'
#!/bin/bash
# Irregular mesh — CPU throughput sweep
set -euo pipefail
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-64}
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
for cfg in test-irregular-3d-*.cfg; do
    echo "--- $cfg ---"
    ./dynearthsol3d "$cfg" 2>&1 | tee "${cfg%.cfg}-cpu.log"
done
echo "Done. Collect: grep 's/step' *-cpu.log"
EOF

cat > "$DIR/run-gpu.sh" << 'EOF'
#!/bin/bash
# Irregular mesh — GPU throughput sweep
set -euo pipefail
for cfg in test-irregular-3d-*.cfg; do
    echo "--- $cfg ---"
    ./dynearthsol3d.gpu "$cfg" 2>&1 | tee "${cfg%.cfg}-gpu.log"
done
echo "Done. Collect: grep 's/step' *-gpu.log"
EOF

chmod +x "$DIR/run-cpu.sh" "$DIR/run-gpu.sh"
echo "  $(ls "$DIR"/*.cfg | wc -l) cfg files, run-cpu.sh + run-gpu.sh"

# ---------------------------------------------------------------------------
# 03-scaling   (regular mesh, CPU core-scaling)
# ---------------------------------------------------------------------------
DIR="$SCRIPT_DIR/03-scaling"
echo "Creating $DIR ..."
mkdir -p "$DIR"

idx=0
for yl in "${SCALING_SIZES[@]}"; do
    idx=$((idx+1)); printf -v n "%02d" $idx
    gen_cfg "$REG_CFG" "$DIR/test-regular-3d-${n}-${yl}.cfg" "$yl" "scale-reg-$yl"
done

copy_exes "$DIR" no   # CPU only

cat > "$DIR/run.sh" << 'EOF'
#!/bin/bash
# Regular mesh — CPU core-scaling test
# Varies OMP_NUM_THREADS across 1 2 4 8 16 32 64 for 4 problem sizes.
set -euo pipefail
for nthreads in 1 2 4 8 16 32 64; do
    export OMP_NUM_THREADS=$nthreads
    echo "=== OMP_NUM_THREADS=$nthreads ==="
    for cfg in test-regular-3d-*.cfg; do
        echo "  $cfg"
        ./dynearthsol3d "$cfg" 2>&1 | tee "${cfg%.cfg}-t${nthreads}.log"
    done
done
echo "Done. Collect: grep 's/step' *-t*.log"
EOF

chmod +x "$DIR/run.sh"
echo "  $(ls "$DIR"/*.cfg | wc -l) cfg files, run.sh"

# ---------------------------------------------------------------------------
echo ""
echo "Setup complete."
echo "  01-regular/   CPU: cd 01-regular  && ./run-cpu.sh"
echo "                GPU: cd 01-regular  && ./run-gpu.sh"
echo "  02-irregular/ CPU: cd 02-irregular && ./run-cpu.sh"
echo "                GPU: cd 02-irregular && ./run-gpu.sh"
echo "  03-scaling/   CPU: cd 03-scaling  && ./run.sh"
