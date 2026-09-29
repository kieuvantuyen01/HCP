#!/usr/bin/env bash
set -e

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo_root"

echo "================================================================="
echo "HCP SOLVER ENVIRONMENT SETUP (PYTHON & RUST)"
echo "================================================================="

# Step 1: System packages (Debian / Ubuntu)
echo "[*] Step 1: Installing system build dependencies..."
if command -v apt-get &> /dev/null; then
    sudo apt-get update -y
    sudo apt-get install -y \
        build-essential \
        gcc \
        g++ \
        make \
        cmake \
        pkg-config \
        libssl-dev \
        libclang-dev \
        clang \
        zlib1g-dev \
        python3 \
        python3-dev \
        python3-pip \
        python3-venv \
        git \
        curl
fi

# Step 2: Pinned runlim used by the final benchmark harness
echo "[*] Step 2: Checking pinned runlim..."
runlim_commit="188f1e07fa233b787589900e0184092b49167706"
runlim_bin="$repo_root/bin/runlim"
runlim_ok=false
if [[ -x "$runlim_bin" && "$("$runlim_bin" --version 2>/dev/null || true)" == "2.0.0rc14" ]]; then
    runlim_ok=true
fi
if [[ "$runlim_ok" != true ]]; then
    echo "[*] Building runlim 2.0.0rc14 from pinned commit $runlim_commit..."
    mkdir -p "$repo_root/bin"
    runlim_build_dir="$(mktemp -d)"
    trap 'rm -rf "$runlim_build_dir"' EXIT
    curl -fsSL \
        "https://github.com/arminbiere/runlim/archive/${runlim_commit}.tar.gz" \
        -o "$runlim_build_dir/runlim.tar.gz"
    tar -xzf "$runlim_build_dir/runlim.tar.gz" -C "$runlim_build_dir"
    runlim_source="$runlim_build_dir/runlim-$runlim_commit"
    (cd "$runlim_source" && ./configure.sh && make)
    install -m 755 "$runlim_source/runlim" "$runlim_bin"
    rm -rf "$runlim_build_dir"
    trap - EXIT
fi
if [[ "$("$runlim_bin" --version)" != "2.0.0rc14" ]]; then
    echo "[!] Unexpected runlim version at $runlim_bin" >&2
    exit 1
fi
echo "[*] runlim ready: $runlim_bin ($("$runlim_bin" --version))"

# Step 3: Rust toolchain (if not already installed)
echo "[*] Step 3: Checking Rust toolchain..."
if ! command -v cargo &> /dev/null; then
    echo "[*] Installing Rust via rustup..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source "$HOME/.cargo/env"
else
    echo "[*] Rust is already installed: $(rustc --version)"
fi

# Build native Rust solver
echo "[*] Building native cegar-fix Rust binary..."
cargo build --release --manifest-path src/cegar-fix/Cargo.toml
echo "[*] Rust binary ready at: src/cegar-fix/target/release/cegar-fix"

# Step 4: Python dependencies
echo "[*] Step 4: Installing Python dependencies..."
python3 -m pip install --upgrade pip --break-system-packages 2>/dev/null || python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt --break-system-packages 2>/dev/null || python3 -m pip install -r requirements.txt

# Step 5: Verification
echo "[*] Step 5: Verifying CaDiCaL SAT Solver in Python..."
python3 -c "from pysat.solvers import Cadical195; s = Cadical195(); s.add_clause([1, 2]); assert s.solve() == True; print('[*] PySAT CaDiCaL 1.9.5: Verified and Working!')"

echo "[*] Step 6: Verifying official 29-graph benchmark..."
python3 scratch/verify_29.py

echo "================================================================="
echo "[*] SUCCESS! Full Python + Rust environment is ready!"
echo "================================================================="
